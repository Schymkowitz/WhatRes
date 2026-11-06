// Clashscore.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// MolProbity-style all-atom clash score.  For every residue we count the
// number of serious steric overlaps (van der Waals overlap >= overlapCut, default
// 0.4 Å) between its atoms and any other atom in the structure.  Hydrogens are
// added geometrically (reusing ZipperScore::AddVirtualHydrogens), so the metric
// is computed identically for every structure regardless of whether the deposited
// model contains hydrogens.
//
// Reuses (from FxAtoms.h): AddVirtualHydrogens, VdwRadius, ElementFromName.
// New code here is only the overlap counting + bonded/H-bond exclusions.
//
// Fidelity note: hydrogens are placed from heavy-atom geometry, NOT optimised
// like Reduce (no Asn/Gln/His flips, no rotatable-H / H-bond-network search).
// So this is "MolProbity-style", not bit-identical to MolProbity's clashscore.
// Heavy-atom and apolar-H clashes are accurate; polar-H clashes are approximate.
// What matters for a resolution-discriminating feature is that it is computed
// *identically* across all structures.
//
// Exclusions (no bond graph needed):
//   - atom pairs within the SAME residue (intra-residue clash is FoldX's
//     torsional_clash term);
//   - the peptide-bond neighbourhood between consecutive residues (1-2 / 1-3
//     pairs around C(i)-N(i+1));
//   - H-bonds: a polar-H ... (N/O/S acceptor) contact is discounted by
//     hbondAllow (default 0.6 Å) before testing the overlap cut.
// -----------------------------------------------------------------------------

#pragma once
#include <string>
#include <vector>

class PDB_File;

struct ClashResult
{
    char   chain      = ' ';
    int    resnum     = 0;
    char   aa1        = '?';
    double clashscore = 0.0;   // # clashing contacts (overlap >= cut) involving
                               // this residue's atoms
    ClashResult() = default;
    ClashResult(char c, int r, char a, double s)
        : chain(c), resnum(r), aa1(a), clashscore(s) {}
};

class Clashscore
{
public:
    // Whole-structure computation: inter-residue clashes need every chain's
    // atoms, so this takes the full PDB_File.  Returns one ClashResult per
    // standard-amino-acid residue, in structure order.
    // includeHydrogens=false (default) → heavy-atom clashscore with the shared
    // covalent bond graph for 1-2/1-3/1-4 exclusion (clean + resolution-
    // discriminating).  includeHydrogens=true adds virtual H, but a faithful
    // all-atom (MolProbity) count needs Reduce-optimised H + Probe radii; the
    // current geometric H over-count normal contacts, so it is OFF by default.
    static std::vector<ClashResult> Compute(const PDB_File& pdb,
                                            double overlapCut = 0.4,
                                            double hbondAllow = 0.6,
                                            bool   includeHydrogens = false);

    // Structure-level clashscore = clashes per 1000 atoms (MolProbity convention).
    // nAtoms = total atoms used (heavy + virtual H).
    static double StructureClashscore(const std::vector<ClashResult>& res,
                                      int nAtoms);
};

// Moved unchanged from Clashscore.cpp (2026-10-07) so other code uses the same radii.
namespace ClashRadii {
// ── MolProbity / probe electron-cloud VdW radii (atomprops.h eRad column) ─────
// Applied ONLY inside this clash metric, so the global ZipperScore::VdwRadius
// used by the XBeta force field is left completely untouched.  probe typing
// distinguishes aromatic ring atoms (Car=1.75, Har=1.05) from aliphatic
// (C=1.70, H=1.22) and polar H (Hpol=1.05); heavy N/O/S are 1.55/1.40/1.80.
// (These are probe's explicit electron-cloud eRad radii, atomprops.h.)
inline bool isAromaticRingAtom(char aa1, const std::string& nm) {
    switch (aa1) {
        case 'F':  // Phe ring C + ring H (CB aliphatic)
            return nm=="CG"||nm=="CD1"||nm=="CD2"||nm=="CE1"||nm=="CE2"||nm=="CZ"
                 ||nm=="HD1"||nm=="HD2"||nm=="HE1"||nm=="HE2"||nm=="HZ";
        case 'Y':  // Tyr ring (CZ bears OH but is aromatic C; HH is O-H, polar)
            return nm=="CG"||nm=="CD1"||nm=="CD2"||nm=="CE1"||nm=="CE2"||nm=="CZ"
                 ||nm=="HD1"||nm=="HD2"||nm=="HE1"||nm=="HE2";
        case 'W':  // Trp indole (NE1 is N; HE1 is N-H, polar)
            return nm=="CG"||nm=="CD1"||nm=="CD2"||nm=="CE2"||nm=="CE3"
                 ||nm=="CZ2"||nm=="CZ3"||nm=="CH2"
                 ||nm=="HD1"||nm=="HE3"||nm=="HZ2"||nm=="HZ3"||nm=="HH2";
        case 'H':  // His imidazole (ND1/NE2 are N; HD1/HE2 are N-H, polar)
            return nm=="CG"||nm=="CD2"||nm=="CE1"||nm=="HD2"||nm=="HE1";
        default: return false;
    }
}

// elem is "C"/"N"/"O"/"S"/"H"; polarH marks an H added as an N/O/S donor proton.
// H radii are probe's explicit electron-cloud eRad column (atomprops.h): the
// nuclear column (1.17/1.00) is NOT what probe's default clashscore uses.
inline float ProbeRadius(char aa1, const std::string& nm, const std::string& elem, bool polarH) {
    if (elem == "H") {
        if (polarH) return 1.05f;                        // Hpol eRad
        if (isAromaticRingAtom(aa1, nm)) return 1.05f;   // Har eRad (aromatic C-H)
        return 1.22f;                                    // aliphatic H eRad
    }
    if (elem == "N") return 1.55f;
    if (elem == "O") return 1.40f;
    if (elem == "S") return 1.80f;
    if (elem == "C") {
        if (nm == "C") return 1.65f;   // backbone carbonyl C (probe O(i)-C(i+1) gaps match only with 1.65)
        return isAromaticRingAtom(aa1, nm) ? 1.75f : 1.70f;  // Car vs C
    }
    return 1.70f;
}
} // namespace ClashRadii
