// SurfacePatcher.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Identifies solvent-exposed residues on a protein chain using a Cβ
// neighbour-count burial model, then clusters exposed residues of the same
// physicochemical type into patches/clusters.
//
// Three patch types are reported:
//   HYDROPHOBIC    residues in `hphob_set` (default: AVILMFWYPCG)
//   POSITIVE       K, R, H
//   NEGATIVE       D, E
//
// All geometric parameters are set as plain public fields before calling
// Analyse(); they map 1:1 to the surface_* entries in Option_class.
//
// Designed to be fast enough for thousands of structures: burial uses a
// SpatialGrid (O(n) build + O(1) per-residue query) and patch adjacency
// uses an O(k²) pairwise check on the small exposed-residue subset.
//
// Usage:
//   SurfacePatcher sp;
//   sp.burial_radius    = O.surface_burial_radius;
//   sp.burial_threshold = O.surface_burial_threshold;
//   // ... etc.
//   SurfaceAnalysis result = sp.Analyse(pdb.peptides[c], pdbname);
//
// -----------------------------------------------------------------------------

#pragma once
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include "Peptide.h"
#include "Residue.h"
#include "Atom.h"
#include "SpatialGrid.h"
#include "vector3D.h"

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
// SurfacePatch: one identified patch or charge cluster
// ─────────────────────────────────────────────────────────────────────────────
struct SurfacePatch
{
    vector<int> indices;        // positions in Peptide::reslist

    // Metrics filled by SurfacePatcher::computeMetrics()
    int    size         = 0;    // number of member residues
    double kd_sum       = 0.0;  // Kyte-Doolittle hydrophobicity sum
    double net_charge   = 0.0;  // net formal charge (K/R=+1, D/E=-1, H=+0.1)
    int    aromatic     = 0;    // count of F, W, Y, H
    double rg           = 0.0;  // radius of gyration of Cβ positions (Å)
    double mean_burial  = 0.0;  // mean Cβ-neighbour count of members
    double comp_entropy = 0.0;  // Shannon entropy of residue-type composition (bits)
    string residue_list;        // e.g. "I42,L45,V48" (one-letter + PDB number)
};

// ─────────────────────────────────────────────────────────────────────────────
// SurfaceAnalysis: complete per-chain result, written as one output row
// ─────────────────────────────────────────────────────────────────────────────
struct SurfaceAnalysis
{
    string pdb_name;
    char   chain            = ' ';
    int    n_total          = 0;    // total residues (ATOM records only)
    int    n_surface        = 0;    // exposed residues (burial < threshold)
    double frac_surface     = 0.0;  // n_surface / n_total
    double mean_burial      = 0.0;  // mean Cβ-neighbour count across all residues
    double surf_net_charge  = 0.0;  // net charge summed over exposed residues only
    double frac_surf_hphob  = 0.0;  // fraction of exposed residues that are hydrophobic

    // Patches sorted by size (largest first)
    vector<SurfacePatch> hydrophobic;
    vector<SurfacePatch> positive;
    vector<SurfacePatch> negative;

    // Convenience accessors: return an empty patch when the vector is empty
    static const SurfacePatch& emptyPatch() { static SurfacePatch ep; return ep; }
    const SurfacePatch& topHydrophobic() const { return hydrophobic.empty() ? emptyPatch() : hydrophobic[0]; }
    const SurfacePatch& topPositive()    const { return positive.empty()    ? emptyPatch() : positive[0];    }
    const SurfacePatch& topNegative()    const { return negative.empty()    ? emptyPatch() : negative[0];    }
};

// ─────────────────────────────────────────────────────────────────────────────
// SurfacePatcher: analysis engine
// ─────────────────────────────────────────────────────────────────────────────
class SurfacePatcher
{
public:
    // --- Geometry (set from Option_class before calling Analyse) -------------
    double burial_radius    = 10.0;         // Cβ-neighbour sphere radius (Å)
    int    burial_threshold = 16;           // exposed if neighbour_count < threshold
    double adj_cutoff       = 8.0;          // Cβ–Cβ distance for patch adjacency (Å)
    int    min_patch_size   = 2;            // minimum residues to report a patch
    double pH               = 7.0;          // used for net-charge calculation
    string hphob_set        = "AVILMFWYPCG"; // one-letter codes classed as hydrophobic

    // Analyse a single chain.  pdb_name is stored in the result for output.
    SurfaceAnalysis Analyse(const Peptide& pep, const string& pdb_name) const;

private:
    // Return the Cβ position for residue r (falls back to Cα for Gly / missing CB).
    Vector3D cbPos(const Residue& r) const;

    // Per-residue classification
    bool   isHydrophobic(char aa) const;
    static bool   isPositive(char aa);
    static bool   isNegative(char aa);
    static double kdValue(char aa);         // Kyte-Doolittle scale
    static double chargeAt7(char aa);       // simplified formal charge at pH 7

    // Build adjacency list (size n) among `candidates` using Cβ-Cβ distance.
    // adj[i] contains local indices j where Cβ-Cβ dist(candidates[i], candidates[j]) < adj_cutoff.
    vector<vector<int>> buildAdjacency(const Peptide& pep,
                                        const vector<int>& candidates) const;

    // BFS connected-components on adjacency list → list of index groups.
    vector<vector<int>> connectedComponents(const vector<vector<int>>& adj, int n) const;

    // Fill all metric fields in a SurfacePatch given its member indices and precomputed burial values.
    SurfacePatch computeMetrics(const Peptide& pep,
                                const vector<int>& indices,
                                const vector<double>& burial) const;
};
