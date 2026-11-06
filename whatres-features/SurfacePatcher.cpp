// SurfacePatcher.cpp
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz

#include "SurfacePatcher.h"
#include <queue>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <functional>

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
// Static physicochemical lookup tables
// ─────────────────────────────────────────────────────────────────────────────

// Kyte-Doolittle hydrophobicity scale (J. Mol. Biol. 157:105-132, 1982).
double SurfacePatcher::kdValue(char aa)
{
    switch (aa)
    {
        case 'A': return  1.8;  case 'R': return -4.5;  case 'N': return -3.5;
        case 'D': return -3.5;  case 'C': return  2.5;  case 'Q': return -3.5;
        case 'E': return -3.5;  case 'G': return -0.4;  case 'H': return -3.2;
        case 'I': return  4.5;  case 'L': return  3.8;  case 'K': return -3.9;
        case 'M': return  1.9;  case 'F': return  2.8;  case 'P': return -1.6;
        case 'S': return -0.8;  case 'T': return -0.7;  case 'W': return -0.9;
        case 'Y': return -1.3;  case 'V': return  4.2;
        default:  return  0.0;
    }
}

// Simplified formal charge at pH 7.
//   K, R fully protonated (+1); D, E fully deprotonated (-1);
//   H ~10% protonated at pH 7 (pKa ~6.0) → +0.1.
double SurfacePatcher::chargeAt7(char aa)
{
    switch (aa)
    {
        case 'K': return  1.0;
        case 'R': return  1.0;
        case 'H': return  0.1;
        case 'D': return -1.0;
        case 'E': return -1.0;
        default:  return  0.0;
    }
}

bool SurfacePatcher::isPositive(char aa) { return aa == 'K' || aa == 'R' || aa == 'H'; }
bool SurfacePatcher::isNegative(char aa) { return aa == 'D' || aa == 'E'; }
bool SurfacePatcher::isHydrophobic(char aa) const
{
    return hphob_set.find(aa) != string::npos;
}

// ─────────────────────────────────────────────────────────────────────────────
// Cβ position (CA fallback for Gly or missing CB)
// ─────────────────────────────────────────────────────────────────────────────
Vector3D SurfacePatcher::cbPos(const Residue& r) const
{
    for (const auto& a : r.atoms)
        if (a.atom_name.find("CB") != string::npos)
            return Vector3D(a.x, a.y, a.z);

    // No CB found: use CA (Gly, or occasionally missing in low-res structures)
    if (r.CA >= 0 && r.CA < (int)r.atoms.size())
        return Vector3D(r.atoms[r.CA].x, r.atoms[r.CA].y, r.atoms[r.CA].z);

    // Last resort: first atom
    return r.atoms.empty() ? Vector3D(0.0, 0.0, 0.0)
                           : Vector3D(r.atoms[0].x, r.atoms[0].y, r.atoms[0].z);
}

// ─────────────────────────────────────────────────────────────────────────────
// Patch adjacency graph
// ─────────────────────────────────────────────────────────────────────────────
vector<vector<int>> SurfacePatcher::buildAdjacency(const Peptide& pep,
                                                    const vector<int>& cands) const
{
    int n = (int)cands.size();
    vector<vector<int>> adj(n);
    const double cut2 = adj_cutoff * adj_cutoff;

    for (int i = 0; i < n; i++)
    {
        Vector3D pi = cbPos(pep.reslist[cands[i]]);
        for (int j = i + 1; j < n; j++)
        {
            Vector3D pj = cbPos(pep.reslist[cands[j]]);
            double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
            if (dx*dx + dy*dy + dz*dz <= cut2)
            {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }
    return adj;
}

// ─────────────────────────────────────────────────────────────────────────────
// BFS connected components
// ─────────────────────────────────────────────────────────────────────────────
vector<vector<int>> SurfacePatcher::connectedComponents(const vector<vector<int>>& adj,
                                                         int n) const
{
    vector<bool>      visited(n, false);
    vector<vector<int>> components;

    for (int start = 0; start < n; start++)
    {
        if (visited[start]) continue;
        vector<int> comp;
        queue<int>  q;
        q.push(start);
        visited[start] = true;
        while (!q.empty())
        {
            int cur = q.front(); q.pop();
            comp.push_back(cur);
            for (int nb : adj[cur])
                if (!visited[nb]) { visited[nb] = true; q.push(nb); }
        }
        components.push_back(comp);
    }
    return components;
}

// ─────────────────────────────────────────────────────────────────────────────
// Patch metric computation
// ─────────────────────────────────────────────────────────────────────────────
SurfacePatch SurfacePatcher::computeMetrics(const Peptide& pep,
                                             const vector<int>& indices,
                                             const vector<double>& burial) const
{
    SurfacePatch p;
    p.indices = indices;
    p.size    = (int)indices.size();

    // Cβ positions for Rg
    vector<Vector3D> positions;
    positions.reserve(p.size);

    // Residue-type counts for Shannon entropy
    map<char, int> typeCount;

    for (int idx : indices)
    {
        const Residue& r = pep.reslist[idx];
        char aa = r.OneLetter();

        p.kd_sum      += kdValue(aa);
        p.net_charge  += chargeAt7(aa);
        p.mean_burial += burial[idx];
        if (aa == 'F' || aa == 'W' || aa == 'Y' || aa == 'H') p.aromatic++;
        typeCount[aa]++;

        positions.push_back(cbPos(r));

        // Build residue list string "A42,I45,..."
        if (!p.residue_list.empty()) p.residue_list += ',';
        p.residue_list += aa;
        p.residue_list += to_string(r.number_in_pdb);
    }

    p.mean_burial /= p.size;

    // Radius of gyration
    double cx = 0, cy = 0, cz = 0;
    for (const auto& v : positions) { cx += v.x; cy += v.y; cz += v.z; }
    cx /= p.size; cy /= p.size; cz /= p.size;
    double sumR2 = 0.0;
    for (const auto& v : positions)
    {
        double dx = v.x - cx, dy = v.y - cy, dz = v.z - cz;
        sumR2 += dx*dx + dy*dy + dz*dz;
    }
    p.rg = sqrt(sumR2 / p.size);

    // Shannon entropy of residue-type composition
    double h = 0.0;
    for (const auto& kv : typeCount)
    {
        double prob = (double)kv.second / p.size;
        h -= prob * log2(prob);
    }
    p.comp_entropy = h;

    return p;
}

// ─────────────────────────────────────────────────────────────────────────────
// Main analysis entry point
// ─────────────────────────────────────────────────────────────────────────────
SurfaceAnalysis SurfacePatcher::Analyse(const Peptide& pep, const string& pdb_name) const
{
    SurfaceAnalysis result;
    result.pdb_name = pdb_name;
    result.chain    = pep.chain;

    // Collect valid (ATOM) residues with at least one heavy atom
    vector<int> valid;
    for (int i = 0; i < (int)pep.reslist.size(); i++)
    {
        const Residue& r = pep.reslist[i];
        if (r.molecule_type == "ATOM" && !r.atoms.empty())
            valid.push_back(i);
    }
    result.n_total = (int)valid.size();
    if (result.n_total == 0) return result;

    // ── Build SpatialGrid of all Cβ positions ─────────────────────────────
    // Cell size = burial_radius ensures the 3×3×3 neighbour-cell search
    // covers the full sphere.  A cell size of burial_radius/2 would halve
    // the per-cell population without requiring extra ring searches.
    SpatialGrid grid(burial_radius / 2.0);
    for (int idx : valid)
        grid.Insert(cbPos(pep.reslist[idx]));

    // ── Per-residue burial depth ──────────────────────────────────────────
    // burial[i] = number of OTHER Cβ atoms within burial_radius of residue i.
    vector<double> burial(pep.reslist.size(), 0.0);
    double totalBurial = 0.0;

    vector<int> surfaceIndices;

    for (int idx : valid)
    {
        Vector3D pos = cbPos(pep.reslist[idx]);
        // CountWithin includes the residue itself, so subtract 1.
        int count = grid.CountWithin(pos, burial_radius) - 1;
        burial[idx] = count;
        totalBurial += count;

        if (count < burial_threshold)
            surfaceIndices.push_back(idx);
    }

    result.mean_burial  = totalBurial / result.n_total;
    result.n_surface    = (int)surfaceIndices.size();
    result.frac_surface = result.n_total > 0 ? (double)result.n_surface / result.n_total : 0.0;

    // ── Global surface charge & hydrophobic fraction ──────────────────────
    int surfHphob = 0;
    for (int idx : surfaceIndices)
    {
        char aa = pep.reslist[idx].OneLetter();
        result.surf_net_charge += chargeAt7(aa);
        if (isHydrophobic(aa)) surfHphob++;
    }
    result.frac_surf_hphob = result.n_surface > 0
                             ? (double)surfHphob / result.n_surface : 0.0;

    // ── Helper: find and score patches of a given type ────────────────────
    auto findPatches = [&](const function<bool(char)>& typeFilter,
                           vector<SurfacePatch>&         patchVec)
    {
        // Candidates: surface residues matching the type filter
        vector<int> cands;
        for (int idx : surfaceIndices)
            if (typeFilter(pep.reslist[idx].OneLetter()))
                cands.push_back(idx);
        if (cands.empty()) return;

        vector<vector<int>> adj   = buildAdjacency(pep, cands);
        vector<vector<int>> comps = connectedComponents(adj, (int)cands.size());

        for (const auto& comp : comps)
        {
            if ((int)comp.size() < min_patch_size) continue;
            // Map local component indices back to peptide indices
            vector<int> pepIndices;
            pepIndices.reserve(comp.size());
            for (int li : comp) pepIndices.push_back(cands[li]);
            patchVec.push_back(computeMetrics(pep, pepIndices, burial));
        }
        // Sort largest first
        sort(patchVec.begin(), patchVec.end(),
             [](const SurfacePatch& a, const SurfacePatch& b){ return a.size > b.size; });
    };

    findPatches([this](char aa){ return isHydrophobic(aa); }, result.hydrophobic);
    findPatches([](char aa){ return isPositive(aa); },        result.positive);
    findPatches([](char aa){ return isNegative(aa); },        result.negative);

    return result;
}
