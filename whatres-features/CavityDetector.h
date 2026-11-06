// CavityDetector.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Detects buried cavities inside a protein structure using a grid-based
// flood-fill approach:
//
//   1. Voxelise: place all heavy atoms onto a 3D grid, mark occupied voxels.
//   2. Expand: dilate occupied voxels by the probe radius (marks protein volume).
//   3. Flood-fill: start from a corner (guaranteed to be solvent), BFS-fill all
//      reachable empty voxels → these become "bulk solvent".
//   4. Whatever remains empty and unreachable is a buried cavity.
//   5. Connected-component labelling on cavity voxels gives individual cavities.
//
// Each cavity is reported with its volume (Å³), centroid, and a list of the
// nearest protein residues lining it.
//
// Usage:
//   CavityDetector cd;
//   cd.grid_spacing  = 1.0;         // Å per voxel (from options)
//   cd.probe_radius  = 1.4;         // water probe
//   cd.min_volume     = 10.0;       // Å³ minimum to report
//   vector<CavityResult> cavities = cd.Detect(pdb);  // all chains combined
//
// The entire PDB_File is passed (not a single Peptide) because cavities
// exist between chains, not within one chain.
// -----------------------------------------------------------------------------

#pragma once
#include <vector>
#include <string>
#include <set>
#include <queue>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include "PDB_File.h"
#include "Peptide.h"
#include "Residue.h"
#include "Atom.h"
#include "vector3D.h"

// ─────────────────────────────────────────────────────────────────────────────
// Result structs
// ─────────────────────────────────────────────────────────────────────────────
struct CavityResult
{
    int    id           = 0;       // cavity index (1-based)
    int    n_voxels     = 0;       // number of cavity voxels
    double volume       = 0.0;     // n_voxels × voxel_volume (Å³)
    double cx           = 0.0;     // centroid x
    double cy           = 0.0;     // centroid y
    double cz           = 0.0;     // centroid z
    std::vector<std::string> lining_residues;   // e.g. "A:V42", "B:L105"
};

struct CavitySummary
{
    int    n_cavities        = 0;       // total cavities found
    double total_volume      = 0.0;     // sum of all cavity volumes (Å³)
    double largest_volume    = 0.0;     // volume of the biggest cavity
    double protein_volume    = 0.0;     // total occupied volume (Å³)
    double void_fraction     = 0.0;     // total_cavity / protein_volume
    std::vector<CavityResult> cavities; // individual cavities, sorted large→small
};

// ─────────────────────────────────────────────────────────────────────────────
// CavityDetector
// ─────────────────────────────────────────────────────────────────────────────
class CavityDetector
{
public:
    // --- Parameters (set from options before calling Detect) -----------------
    double grid_spacing = 1.0;      // Å per voxel edge
    double probe_radius = 1.4;      // solvent probe radius (Å)
    double min_volume   = 10.0;     // minimum cavity volume to report (Å³)
    double lining_dist  = 3.5;      // max distance from cavity voxel to atom to be "lining"

    // Detect cavities in the entire PDB structure (all chains).
    CavitySummary Detect(const PDB_File& pdb) const
    {
        CavitySummary summary;

        // ── Collect all heavy atoms with chain+residue info ──────────────
        struct AtomInfo {
            double x, y, z;
            char   chain;
            char   aa1;
            int    resnum;
        };
        std::vector<AtomInfo> atoms;
        atoms.reserve(5000);

        for (const auto& pep : pdb.peptides)
        {
            for (const auto& res : pep.reslist)
            {
                if (res.molecule_type != "ATOM") continue;
                char aa1 = res.OneLetter();
                for (const auto& a : res.atoms)
                {
                    // Skip hydrogens
                    if (!a.atom_name.empty() && a.atom_name[0] == 'H') continue;
                    if (a.atom_name.size() > 1 && a.atom_name[0] == ' ' && a.atom_name[1] == 'H') continue;
                    atoms.push_back({a.x, a.y, a.z, res.chain, aa1, res.number_in_pdb});
                }
            }
        }

        if (atoms.empty()) return summary;

        // ── Bounding box with padding ────────────────────────────────────
        double pad = probe_radius + grid_spacing * 3;   // extra margin
        double xmin = atoms[0].x, xmax = atoms[0].x;
        double ymin = atoms[0].y, ymax = atoms[0].y;
        double zmin = atoms[0].z, zmax = atoms[0].z;
        for (const auto& a : atoms)
        {
            if (a.x < xmin) xmin = a.x; if (a.x > xmax) xmax = a.x;
            if (a.y < ymin) ymin = a.y; if (a.y > ymax) ymax = a.y;
            if (a.z < zmin) zmin = a.z; if (a.z > zmax) zmax = a.z;
        }
        xmin -= pad; ymin -= pad; zmin -= pad;
        xmax += pad; ymax += pad; zmax += pad;

        int nx = (int)std::ceil((xmax - xmin) / grid_spacing) + 1;
        int ny = (int)std::ceil((ymax - ymin) / grid_spacing) + 1;
        int nz = (int)std::ceil((zmax - zmin) / grid_spacing) + 1;

        // ── Voxel states: 0=empty, 1=protein, 2=bulk_solvent, 3=cavity ──
        // Use flat array for speed
        long long totalVoxels = (long long)nx * ny * nz;
        if (totalVoxels > 200000000LL)
        {
            // Safety: for very large structures at fine grid, refuse gracefully
            // (200M voxels ~ 200 MB)
            return summary;
        }
        std::vector<unsigned char> grid(totalVoxels, 0);

        auto idx = [&](int ix, int iy, int iz) -> long long {
            return (long long)ix * ny * nz + (long long)iy * nz + iz;
        };

        // ── Mark protein voxels (atom + probe expansion) ─────────────────
        // For each atom, mark all voxels within probe_radius as protein.
        double expand = probe_radius;
        int nExpand = (int)std::ceil(expand / grid_spacing) + 1;

        for (const auto& a : atoms)
        {
            int cix = (int)std::floor((a.x - xmin) / grid_spacing);
            int ciy = (int)std::floor((a.y - ymin) / grid_spacing);
            int ciz = (int)std::floor((a.z - zmin) / grid_spacing);

            for (int dx = -nExpand; dx <= nExpand; dx++)
            for (int dy = -nExpand; dy <= nExpand; dy++)
            for (int dz = -nExpand; dz <= nExpand; dz++)
            {
                int ix = cix + dx, iy = ciy + dy, iz = ciz + dz;
                if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
                    continue;
                // Check actual distance from voxel centre to atom
                double vx = xmin + ix * grid_spacing;
                double vy = ymin + iy * grid_spacing;
                double vz = zmin + iz * grid_spacing;
                double d2 = (vx-a.x)*(vx-a.x) + (vy-a.y)*(vy-a.y) + (vz-a.z)*(vz-a.z);
                if (d2 <= expand * expand)
                    grid[idx(ix, iy, iz)] = 1;  // protein
            }
        }

        // Count protein voxels
        long long nProtein = 0;
        for (auto v : grid) if (v == 1) nProtein++;
        double voxVol = grid_spacing * grid_spacing * grid_spacing;
        summary.protein_volume = nProtein * voxVol;

        // ── Flood-fill from corners → mark bulk solvent ──────────────────
        // All border voxels that are empty are seeds.
        std::queue<long long> q;
        auto tryEnqueue = [&](int ix, int iy, int iz) {
            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
                return;
            long long i = idx(ix, iy, iz);
            if (grid[i] == 0)
            {
                grid[i] = 2;  // bulk solvent
                q.push(i);
            }
        };

        // Seed all border faces
        for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
            {
                tryEnqueue(ix, iy, 0);
                tryEnqueue(ix, iy, nz - 1);
            }
        for (int ix = 0; ix < nx; ix++)
            for (int iz = 0; iz < nz; iz++)
            {
                tryEnqueue(ix, 0, iz);
                tryEnqueue(ix, ny - 1, iz);
            }
        for (int iy = 0; iy < ny; iy++)
            for (int iz = 0; iz < nz; iz++)
            {
                tryEnqueue(0, iy, iz);
                tryEnqueue(nx - 1, iy, iz);
            }

        // BFS flood fill (6-connected)
        while (!q.empty())
        {
            long long ci = q.front(); q.pop();
            // Decode index
            int iz = ci % nz;
            int iy = (ci / nz) % ny;
            int ix = ci / ((long long)ny * nz);

            tryEnqueue(ix - 1, iy, iz);
            tryEnqueue(ix + 1, iy, iz);
            tryEnqueue(ix, iy - 1, iz);
            tryEnqueue(ix, iy + 1, iz);
            tryEnqueue(ix, iy, iz - 1);
            tryEnqueue(ix, iy, iz + 1);
        }

        // ── Remaining empty voxels (state 0) are cavity voxels ───────────
        // Connected-component labelling
        int cavityId = 0;
        std::vector<int> label(totalVoxels, 0);  // 0 = unlabelled

        for (long long ci = 0; ci < totalVoxels; ci++)
        {
            if (grid[ci] != 0) continue;    // only empty (cavity) voxels
            if (label[ci] != 0) continue;   // already labelled

            cavityId++;
            std::queue<long long> bfs;
            bfs.push(ci);
            label[ci] = cavityId;

            while (!bfs.empty())
            {
                long long cur = bfs.front(); bfs.pop();
                int iz2 = cur % nz;
                int iy2 = (cur / nz) % ny;
                int ix2 = cur / ((long long)ny * nz);

                auto tryLabel = [&](int ax, int ay, int az) {
                    if (ax < 0 || ax >= nx || ay < 0 || ay >= ny || az < 0 || az >= nz)
                        return;
                    long long ni = idx(ax, ay, az);
                    if (grid[ni] == 0 && label[ni] == 0)
                    {
                        label[ni] = cavityId;
                        bfs.push(ni);
                    }
                };

                tryLabel(ix2 - 1, iy2, iz2);
                tryLabel(ix2 + 1, iy2, iz2);
                tryLabel(ix2, iy2 - 1, iz2);
                tryLabel(ix2, iy2 + 1, iz2);
                tryLabel(ix2, iy2, iz2 - 1);
                tryLabel(ix2, iy2, iz2 + 1);
            }
        }

        // ── Collect cavity statistics ────────────────────────────────────
        // Count voxels + centroid per cavity
        struct CavAcc { int count = 0; double sx = 0, sy = 0, sz = 0; };
        std::unordered_map<int, CavAcc> cavMap;

        for (long long ci = 0; ci < totalVoxels; ci++)
        {
            if (label[ci] == 0) continue;
            int iz2 = ci % nz;
            int iy2 = (ci / nz) % ny;
            int ix2 = ci / ((long long)ny * nz);
            double vx = xmin + ix2 * grid_spacing;
            double vy = ymin + iy2 * grid_spacing;
            double vz = zmin + iz2 * grid_spacing;
            auto& ca = cavMap[label[ci]];
            ca.count++;
            ca.sx += vx; ca.sy += vy; ca.sz += vz;
        }

        // ── Build CavityResult for each cavity above min_volume ──────────
        for (auto& kv : cavMap)
        {
            double vol = kv.second.count * voxVol;
            if (vol < min_volume) continue;

            CavityResult cr;
            cr.id       = kv.first;
            cr.n_voxels = kv.second.count;
            cr.volume   = vol;
            cr.cx       = kv.second.sx / kv.second.count;
            cr.cy       = kv.second.sy / kv.second.count;
            cr.cz       = kv.second.sz / kv.second.count;

            // Find lining residues: residues with any atom within lining_dist
            // of the cavity centroid.  For large cavities we sample a few
            // representative voxels instead of all, but centroid is sufficient
            // for the residue list.
            std::set<std::string> seen;
            for (const auto& a : atoms)
            {
                double d2 = (a.x - cr.cx)*(a.x - cr.cx) +
                            (a.y - cr.cy)*(a.y - cr.cy) +
                            (a.z - cr.cz)*(a.z - cr.cz);
                // Use a generous radius: half the cavity "diameter" + lining_dist
                double cavRadius = std::cbrt(vol * 3.0 / (4.0 * 3.14159)) + lining_dist;
                if (d2 <= cavRadius * cavRadius)
                {
                    std::string key;
                    key += a.chain;
                    key += ':';
                    key += a.aa1;
                    key += std::to_string(a.resnum);
                    if (seen.insert(key).second)
                        cr.lining_residues.push_back(key);
                }
            }

            summary.cavities.push_back(cr);
        }

        // Sort by volume, largest first
        std::sort(summary.cavities.begin(), summary.cavities.end(),
                  [](const CavityResult& a, const CavityResult& b) {
                      return a.volume > b.volume;
                  });

        // Renumber 1-based
        for (int i = 0; i < (int)summary.cavities.size(); i++)
            summary.cavities[i].id = i + 1;

        summary.n_cavities = (int)summary.cavities.size();
        summary.total_volume = 0;
        for (const auto& c : summary.cavities)
            summary.total_volume += c.volume;
        summary.largest_volume = summary.cavities.empty() ? 0.0 : summary.cavities[0].volume;
        summary.void_fraction  = summary.protein_volume > 0
                                 ? summary.total_volume / summary.protein_volume : 0.0;

        return summary;
    }

};
