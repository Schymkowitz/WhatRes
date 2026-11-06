// InterfaceVoidAnalyser.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Quantifies empty space at the interface between two chains in a complex.
//
// The problem: in protein-protein complexes (especially antibody-antigen),
// the interface is not fully enclosed, so flood-fill cavity detection marks
// it as bulk solvent.  But comparing how tightly two chains pack (how
// much "void" exists between them) is critical for antibody engineering.
//
// Approach:
//   1. Identify interface residues: residues on chain A within contact_dist
//      of any atom on chain B, and vice versa.
//   2. Define an interface bounding box from these residues (with padding).
//   3. Voxelise this box:  mark voxels occupied by either chain.
//   4. For each empty voxel, check if it sits "between" the two chains:
//      it must have at least one atom from chain A within shell_dist AND
//      at least one atom from chain B within shell_dist.
//   5. These "interface void" voxels are the empty space at the interface.
//
// Metrics reported:
//   - interface_area_approx   (number of interface residues × ~area per residue)
//   - interface_void_volume   (Å³ of empty voxels between the two chains)
//   - interface_packing_density (1 - void_fraction in the interface box)
//   - interface_void_index    (void_volume / n_interface_residues, normalised)
//   - interface_residues_A, interface_residues_B  (residue lists)
//
// For Rob's use case (comparing many variants of the same antibody-antigen
// complex), the void_index gives a single number for how well-packed the
// interface is, directly comparable across mutants.
//
// Usage:
//   InterfaceVoidAnalyser iva;
//   iva.grid_spacing = 0.5;
//   InterfaceVoidResult r = iva.Analyse(pdb, 'H', 'A');  // heavy chain vs antigen
//
// -----------------------------------------------------------------------------

#pragma once
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "PDB_File.h"
#include "Peptide.h"
#include "Residue.h"
#include "Atom.h"
#include "vector3D.h"

// ─────────────────────────────────────────────────────────────────────────────
// Result
// ─────────────────────────────────────────────────────────────────────────────
struct InterfaceVoidResult
{
    char   chainA           = ' ';
    char   chainB           = ' ';
    int    n_interface_A    = 0;    // interface residues on chain A
    int    n_interface_B    = 0;    // interface residues on chain B
    double interface_void_volume = 0.0;   // Å³
    double interface_total_volume = 0.0;  // total volume of interface box (Å³)
    double interface_occupied     = 0.0;  // occupied volume (Å³)
    double packing_density       = 0.0;   // occupied / total in interface zone
    double void_index            = 0.0;   // void_volume / n_interface_residues
    std::string residues_A;        // e.g. "Y33,R50,D54"
    std::string residues_B;        // e.g. "K417,E484,N501"
};

// ─────────────────────────────────────────────────────────────────────────────
// InterfaceVoidAnalyser
// ─────────────────────────────────────────────────────────────────────────────
class InterfaceVoidAnalyser
{
public:
    // --- Parameters (set from options) ───────────────────────────────────────
    double grid_spacing  = 0.5;    // Å per voxel (finer than cavity detection)
    double contact_dist  = 5.0;    // Å cutoff to define interface residues
    double shell_dist    = 4.0;    // Å; void voxel must have both chains within this
    double probe_radius  = 1.4;    // solvent probe for marking occupied voxels
    double box_padding   = 3.0;    // Å padding around interface residues

    // Analyse the interface between two specific chains.
    // Returns empty result if either chain is missing.
    InterfaceVoidResult Analyse(const PDB_File& pdb, char chainA, char chainB) const
    {
        InterfaceVoidResult result;
        result.chainA = chainA;
        result.chainB = chainB;

        // ── Find the two chains ──────────────────────────────────────────
        const Peptide* pepA = nullptr;
        const Peptide* pepB = nullptr;
        for (const auto& p : pdb.peptides)
        {
            if (p.chain == chainA) pepA = &p;
            if (p.chain == chainB) pepB = &p;
        }
        if (!pepA || !pepB) return result;

        // ── Collect atoms per chain ──────────────────────────────────────
        struct AInfo {
            double x, y, z;
            int    resIdx;     // index in peptide's reslist
            int    chainSide;  // 0 = A, 1 = B
        };
        std::vector<AInfo> atomsA, atomsB;

        auto collectAtoms = [](const Peptide& pep, int side, std::vector<AInfo>& out) {
            for (int ri = 0; ri < (int)pep.reslist.size(); ri++)
            {
                const Residue& r = pep.reslist[ri];
                if (r.molecule_type != "ATOM") continue;
                for (const auto& a : r.atoms)
                {
                    if (!a.atom_name.empty() && a.atom_name[0] == 'H') continue;
                    if (a.atom_name.size() > 1 && a.atom_name[0] == ' ' && a.atom_name[1] == 'H') continue;
                    out.push_back({a.x, a.y, a.z, ri, side});
                }
            }
        };
        collectAtoms(*pepA, 0, atomsA);
        collectAtoms(*pepB, 1, atomsB);

        if (atomsA.empty() || atomsB.empty()) return result;

        // ── Identify interface residues ──────────────────────────────────
        // A residue on chain A is "interface" if any of its atoms is within
        // contact_dist of any atom on chain B, and vice versa.
        double cd2 = contact_dist * contact_dist;
        std::set<int> ifaceA, ifaceB;

        // Use a simple spatial hash on chain B for speed
        double cellSize = contact_dist;
        struct CellKey {
            int x, y, z;
            bool operator==(const CellKey& o) const { return x==o.x && y==o.y && z==o.z; }
        };
        struct CellHash {
            size_t operator()(const CellKey& k) const {
                return std::hash<long long>()(
                    ((long long)k.x << 42) ^ ((long long)k.y << 21) ^ (long long)k.z);
            }
        };
        std::unordered_map<CellKey, std::vector<int>, CellHash> gridB;
        for (int i = 0; i < (int)atomsB.size(); i++)
        {
            CellKey k{(int)std::floor(atomsB[i].x / cellSize),
                      (int)std::floor(atomsB[i].y / cellSize),
                      (int)std::floor(atomsB[i].z / cellSize)};
            gridB[k].push_back(i);
        }

        for (const auto& aa : atomsA)
        {
            int cix = (int)std::floor(aa.x / cellSize);
            int ciy = (int)std::floor(aa.y / cellSize);
            int ciz = (int)std::floor(aa.z / cellSize);
            bool found = false;
            for (int dx = -1; dx <= 1 && !found; dx++)
            for (int dy = -1; dy <= 1 && !found; dy++)
            for (int dz = -1; dz <= 1 && !found; dz++)
            {
                auto it = gridB.find({cix+dx, ciy+dy, ciz+dz});
                if (it == gridB.end()) continue;
                for (int bi : it->second)
                {
                    double ddx = aa.x - atomsB[bi].x;
                    double ddy = aa.y - atomsB[bi].y;
                    double ddz = aa.z - atomsB[bi].z;
                    if (ddx*ddx + ddy*ddy + ddz*ddz <= cd2)
                    {
                        ifaceA.insert(aa.resIdx);
                        ifaceB.insert(atomsB[bi].resIdx);
                        found = true;
                        break;
                    }
                }
            }
        }

        result.n_interface_A = (int)ifaceA.size();
        result.n_interface_B = (int)ifaceB.size();

        // Build residue list strings
        auto buildResList = [](const Peptide& pep, const std::set<int>& indices) -> std::string {
            std::string s;
            for (int idx : indices)
            {
                if (!s.empty()) s += ',';
                s += pep.reslist[idx].OneLetter();
                s += std::to_string(pep.reslist[idx].number_in_pdb);
            }
            return s;
        };
        result.residues_A = buildResList(*pepA, ifaceA);
        result.residues_B = buildResList(*pepB, ifaceB);

        if (ifaceA.empty() || ifaceB.empty()) return result;

        // ── Interface bounding box ───────────────────────────────────────
        // Collect all atoms from interface residues on both chains
        std::vector<AInfo> ifAtoms;
        for (const auto& a : atomsA)
            if (ifaceA.count(a.resIdx)) ifAtoms.push_back(a);
        for (const auto& a : atomsB)
            if (ifaceB.count(a.resIdx)) ifAtoms.push_back(a);

        double xmin = ifAtoms[0].x, xmax = ifAtoms[0].x;
        double ymin = ifAtoms[0].y, ymax = ifAtoms[0].y;
        double zmin = ifAtoms[0].z, zmax = ifAtoms[0].z;
        for (const auto& a : ifAtoms)
        {
            if (a.x < xmin) xmin = a.x; if (a.x > xmax) xmax = a.x;
            if (a.y < ymin) ymin = a.y; if (a.y > ymax) ymax = a.y;
            if (a.z < zmin) zmin = a.z; if (a.z > zmax) zmax = a.z;
        }
        xmin -= box_padding; ymin -= box_padding; zmin -= box_padding;
        xmax += box_padding; ymax += box_padding; zmax += box_padding;

        int nx = (int)std::ceil((xmax - xmin) / grid_spacing) + 1;
        int ny = (int)std::ceil((ymax - ymin) / grid_spacing) + 1;
        int nz = (int)std::ceil((zmax - zmin) / grid_spacing) + 1;

        long long totalVoxels = (long long)nx * ny * nz;
        if (totalVoxels > 100000000LL) return result;  // safety cap

        // ── Voxelise: mark occupied voxels ───────────────────────────────
        // 0 = empty, 1 = occupied by chain A, 2 = occupied by chain B, 3 = both
        std::vector<unsigned char> vox(totalVoxels, 0);

        auto vidx = [&](int ix, int iy, int iz) -> long long {
            return (long long)ix * ny * nz + (long long)iy * nz + iz;
        };

        int nExp = (int)std::ceil(probe_radius / grid_spacing) + 1;
        double pr2 = probe_radius * probe_radius;

        auto markAtoms = [&](const std::vector<AInfo>& atoms_list, unsigned char bit) {
            for (const auto& a : atoms_list)
            {
                int cix = (int)std::floor((a.x - xmin) / grid_spacing);
                int ciy = (int)std::floor((a.y - ymin) / grid_spacing);
                int ciz = (int)std::floor((a.z - zmin) / grid_spacing);
                for (int dx = -nExp; dx <= nExp; dx++)
                for (int dy = -nExp; dy <= nExp; dy++)
                for (int dz = -nExp; dz <= nExp; dz++)
                {
                    int ix = cix+dx, iy = ciy+dy, iz = ciz+dz;
                    if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
                        continue;
                    double vx = xmin + ix * grid_spacing;
                    double vy = ymin + iy * grid_spacing;
                    double vz = zmin + iz * grid_spacing;
                    double d2 = (vx-a.x)*(vx-a.x)+(vy-a.y)*(vy-a.y)+(vz-a.z)*(vz-a.z);
                    if (d2 <= pr2)
                        vox[vidx(ix, iy, iz)] |= bit;
                }
            }
        };

        // Mark ALL atoms from both chains (not just interface), because
        // non-interface atoms can still occupy space in the interface box.
        markAtoms(atomsA, 1);  // bit 0 = chain A
        markAtoms(atomsB, 2);  // bit 1 = chain B

        // ── Count interface void voxels ──────────────────────────────────
        // A void voxel is: empty (vox==0) AND has at least one chain A atom
        // within shell_dist AND at least one chain B atom within shell_dist.
        //
        // For speed, build spatial hashes of both atom sets.
        double sd2 = shell_dist * shell_dist;

        // Simple approach: for each empty voxel, check nearest atom on each chain.
        // This is O(V × A) in the worst case but the interface box is small and
        // we use the spatial hash.
        struct ShellGrid {
            double cell;
            std::unordered_map<long long, std::vector<Vector3D>> bins;
            void insert(double x, double y, double z) {
                int ix = (int)std::floor(x/cell), iy = (int)std::floor(y/cell), iz = (int)std::floor(z/cell);
                bins[((long long)ix<<42)^((long long)iy<<21)^(long long)iz].push_back(Vector3D(x,y,z));
            }
            bool anyWithin(double px, double py, double pz, double dist2) const {
                int ix=(int)std::floor(px/cell), iy=(int)std::floor(py/cell), iz=(int)std::floor(pz/cell);
                for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++) {
                    auto it=bins.find(((long long)(ix+dx)<<42)^((long long)(iy+dy)<<21)^(long long)(iz+dz));
                    if(it==bins.end()) continue;
                    for(const auto& q:it->second) {
                        double ddx=px-q.x, ddy=py-q.y, ddz=pz-q.z;
                        if(ddx*ddx+ddy*ddy+ddz*ddz<=dist2) return true;
                    }
                }
                return false;
            }
        };

        ShellGrid sgA{shell_dist}, sgB{shell_dist};
        for (const auto& a : atomsA) sgA.insert(a.x, a.y, a.z);
        for (const auto& a : atomsB) sgB.insert(a.x, a.y, a.z);

        long long nEmpty = 0, nVoid = 0, nOccupied = 0;
        for (int ix = 0; ix < nx; ix++)
        for (int iy = 0; iy < ny; iy++)
        for (int iz = 0; iz < nz; iz++)
        {
            long long vi = vidx(ix, iy, iz);
            if (vox[vi] != 0)
            {
                nOccupied++;
                continue;
            }
            nEmpty++;

            double vx = xmin + ix * grid_spacing;
            double vy = ymin + iy * grid_spacing;
            double vz = zmin + iz * grid_spacing;

            if (sgA.anyWithin(vx, vy, vz, sd2) && sgB.anyWithin(vx, vy, vz, sd2))
                nVoid++;
        }

        double voxVol = grid_spacing * grid_spacing * grid_spacing;
        result.interface_void_volume  = nVoid * voxVol;
        result.interface_total_volume = totalVoxels * voxVol;
        result.interface_occupied     = nOccupied * voxVol;
        result.packing_density        = totalVoxels > 0
                                        ? (double)nOccupied / totalVoxels : 0.0;
        int nIfaceRes = result.n_interface_A + result.n_interface_B;
        result.void_index = nIfaceRes > 0
                            ? result.interface_void_volume / nIfaceRes : 0.0;

        return result;
    }

    // Analyse ALL chain pairs in a PDB (for multi-chain complexes).
    std::vector<InterfaceVoidResult> AnalyseAll(const PDB_File& pdb) const
    {
        std::vector<InterfaceVoidResult> results;
        int nChains = (int)pdb.peptides.size();
        for (int i = 0; i < nChains; i++)
        for (int j = i + 1; j < nChains; j++)
        {
            char cA = pdb.peptides[i].chain;
            char cB = pdb.peptides[j].chain;
            InterfaceVoidResult r = Analyse(pdb, cA, cB);
            if (r.n_interface_A > 0 && r.n_interface_B > 0)
                results.push_back(r);
        }
        return results;
    }
};
