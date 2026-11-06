// FxAtoms.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Atom-level helpers for the structure-quality commands, taken unchanged from
// ZipperScorer.h of the Switch Lab structure toolkit:
//
//   VdwRadius, ElementFromName      Bondi radii and element from a PDB atom name
//   PlaceAtom, Dist2                internal-coordinate placement, squared distance
//   FxAtomParam, FxAtom             an atom with its FoldX solvation parameters
//   BuildResidueAdjacency           covalent bond graph of one residue (Clashscore)
//   AddVirtualHydrogens             geometric hydrogen placement (Clashscore)
//   FxLookup                        per-atom volume, Occ and Occmax of the FoldX 2014
//                                   parameter set, from the compiled library
//                                   libatomparams.a (atomparams.h)
//   ComputeOccupancy                FoldX occupancy, the Gaussian-weighted volume of
//                                   neighbouring atoms (NativeQCInputs burial)
//
// The atom parameters and the occupancy method are those of FoldX:
//   Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L.
//   The FoldX web server: an online force field.
//   Nucleic Acids Res. 2005;33:W382-W388.
#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>

#include "vector3D.h"
#include "atomparams.h"

namespace ZipperScore {

// ── Van der Waals radii (Å) for heavy atoms ────────────────────────────
// Bondi radii, standard for clash evaluation
inline double VdwRadius(const std::string& elem) {
    if (elem == "C") return 1.70;
    if (elem == "N") return 1.55;
    if (elem == "O") return 1.52;
    if (elem == "S") return 1.80;
    return 1.70;  // default: carbon-like
}

// Element from PDB atom name
inline std::string ElementFromName(const std::string& name) {
    if (name == "N" || name == "NE" || name == "NH1" || name == "NH2" ||
        name == "ND1" || name == "ND2" || name == "NE1" || name == "NE2" ||
        name == "NZ")
        return "N";
    if (name == "O" || name == "OG" || name == "OG1" || name == "OD1" ||
        name == "OD2" || name == "OE1" || name == "OE2" || name == "OH")
        return "O";
    if (name == "S" || name == "SG" || name == "SD")
        return "S";
    return "C";  // CA, CB, CG, CD, CE, CZ, CH2, ...
}

// ── Sidechain atom placement ───────────────────────────────────────────

/// Place atom D given atoms A-B-C, bond length C-D, angle B-C-D, dihedral A-B-C-D
/// This is the fundamental operation for building sidechains from dihedrals.
inline Vector3D PlaceAtom(const Vector3D& A, const Vector3D& B, const Vector3D& C,
                          double bondLen, double angleDeg, double dihedralDeg)
{
    // Convert to radians
    double angle = angleDeg * M_PI / 180.0;
    double dihed = dihedralDeg * M_PI / 180.0;

    // Local coordinate system at C looking back along B-C
    Vector3D bc = C - B;
    double bc_len = std::sqrt(bc.x*bc.x + bc.y*bc.y + bc.z*bc.z);
    if (bc_len < 1e-8) return C;  // degenerate
    Vector3D n = bc;
    n.x /= bc_len; n.y /= bc_len; n.z /= bc_len;

    Vector3D ab = B - A;
    // d = ab × n  (perpendicular in the A-B-C plane)
    Vector3D d = ab % n;
    double d_len = std::sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
    if (d_len < 1e-8) return C;
    d.x /= d_len; d.y /= d_len; d.z /= d_len;

    // m = n × d  (completes right-handed frame)
    Vector3D m = n % d;

    // Place D in local frame then transform to global
    double sa = std::sin(angle);
    double ca_val = std::cos(angle);
    double sd = std::sin(dihed);
    double cd = std::cos(dihed);

    Vector3D D;
    D.x = C.x + bondLen * (-n.x * ca_val + m.x * sa * cd + d.x * sa * sd);
    D.y = C.y + bondLen * (-n.y * ca_val + m.y * sa * cd + d.y * sa * sd);
    D.z = C.z + bondLen * (-n.z * ca_val + m.z * sa * cd + d.z * sa * sd);
    return D;
}


/// Distance squared between two vectors
inline double Dist2(const Vector3D& a, const Vector3D& b) {
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

/// FoldX-style atom parameters (from parameters_solvation.cpp)
struct FxAtomParam {
    double volume;     // atomic volume (Å³)
    double Occ;        // min occupancy (unfolded state)
    double Occmax;     // max occupancy (fully buried)
    double vdw;        // VdW energy (unscaled, kcal/mol)
    double enerG;      // solvation energy (neg = hydrophobic, pos = polar)
    double radius;     // VdW radius for contacts (Å)
    double min_radius; // VdW radius for clashes (Å)
    int    polar;      // 1 = polar atom
    int    nDonor;     // number of H-bond donor protons (0 for non-donors)
    int    nAcceptor;  // number of H-bond acceptor lone pairs (0 for non-acceptors)
};

/// Extended atom with FoldX parameters attached
struct FxAtom {
    std::string name;
    Vector3D    pos;
    FxAtomParam param;
    double      charge = 0.0;   // formal charge (e.g. -0.5 for Asp OD1, +1.0 for Lys NZ)
    int         resSeqIdx = -1; // residue sequential index within chain (-1 = unset)
    // Covalently-bonded heavy atom defining this polar atom's H-bond geometry
    // (e.g. C for backbone O, CD for Asn OD1/ND2). Used for angle-dependent
    // H-bond scoring. Populated once, from residue topology, in FinalizeAtoms.
    Vector3D    hbBase;
    bool        hasHbBase = false;
    // Aromatic ring-plane normal (unit vector), stamped onto Phe/Tyr/Trp/His
    // ring atoms in AddAromaticCentroids. Enables an anisotropic (oblate) vdW
    // radius for ring atoms: the effective clash radius is smaller perpendicular
    // to the ring face than in-plane, a shape measured directly from
    // high-resolution structures. Unset (isAromRing=false) for all other atoms.
    Vector3D    ringNormal;
    bool        isAromRing = false;
};


// Build the covalent bond-adjacency list for ONE residue's atoms (heavy + any
// virtual H already added): the definition of connectivity Clashscore uses
// for 1-2 / 1-3 / 1-4 exclusion.
//   adj[i] = indices of atoms covalently bonded to atom i.
// Returns false (no adjacency built) for unknown residue types or if CA is
// missing.  GLY (backbone only) and PRO (with the CD-N ring bond) are handled.
inline bool BuildResidueAdjacency(const std::vector<FxAtom>& atoms, char aa,
                                  std::vector<std::vector<int>>& adj) {
    using BondList = std::vector<std::pair<std::string, std::string>>;
    static const BondList bbBonds = {{"N","CA"},{"CA","C"},{"C","O"},{"CA","CB"}};

    BondList scBonds;
    switch (aa) {
        case 'A': break; // no SC bonds beyond CB
        case 'G': break; // backbone only (no CB)
        case 'P': scBonds = {{"CB","CG"},{"CG","CD"},{"CD","N"}}; break; // ring
        case 'V': scBonds = {{"CB","CG1"},{"CB","CG2"}}; break;
        case 'I': scBonds = {{"CB","CG1"},{"CB","CG2"},{"CG1","CD1"}}; break;
        case 'L': scBonds = {{"CB","CG"},{"CG","CD1"},{"CG","CD2"}}; break;
        case 'S': scBonds = {{"CB","OG"}}; break;
        case 'T': scBonds = {{"CB","OG1"},{"CB","CG2"}}; break;
        case 'C': scBonds = {{"CB","SG"}}; break;
        case 'M': scBonds = {{"CB","CG"},{"CG","SD"},{"SD","CE"}}; break;
        case 'F': scBonds = {{"CB","CG"},{"CG","CD1"},{"CG","CD2"},
                             {"CD1","CE1"},{"CD2","CE2"},{"CE1","CZ"},{"CE2","CZ"}}; break;
        case 'Y': scBonds = {{"CB","CG"},{"CG","CD1"},{"CG","CD2"},
                             {"CD1","CE1"},{"CD2","CE2"},{"CE1","CZ"},{"CE2","CZ"},
                             {"CZ","OH"}}; break;
        case 'W': scBonds = {{"CB","CG"},{"CG","CD1"},{"CG","CD2"},
                             {"CD1","NE1"},{"CD2","CE2"},{"NE1","CE2"},
                             {"CD2","CE3"},{"CE2","CZ2"},{"CE3","CZ3"},
                             {"CZ2","CH2"},{"CZ3","CH2"}}; break;
        case 'H': scBonds = {{"CB","CG"},{"CG","ND1"},{"CG","CD2"},
                             {"ND1","CE1"},{"CD2","NE2"},{"CE1","NE2"}}; break;
        case 'D': scBonds = {{"CB","CG"},{"CG","OD1"},{"CG","OD2"}}; break;
        case 'N': scBonds = {{"CB","CG"},{"CG","OD1"},{"CG","ND2"}}; break;
        case 'E': scBonds = {{"CB","CG"},{"CG","CD"},{"CD","OE1"},{"CD","OE2"}}; break;
        case 'Q': scBonds = {{"CB","CG"},{"CG","CD"},{"CD","OE1"},{"CD","NE2"}}; break;
        case 'K': scBonds = {{"CB","CG"},{"CG","CD"},{"CD","CE"},{"CE","NZ"}}; break;
        case 'R': scBonds = {{"CB","CG"},{"CG","CD"},{"CD","NE"},{"NE","CZ"},
                             {"CZ","NH1"},{"CZ","NH2"}}; break;
        default: return false;
    }

    // Map atom name → index for the atoms we actually have
    std::map<std::string, int> nameIdx;
    for (int i = 0; i < (int)atoms.size(); ++i) {
        if (atoms[i].name == "X_PI") continue; // skip virtual centroids
        nameIdx[atoms[i].name] = i;
    }
    if (nameIdx.find("CA") == nameIdx.end()) return false; // minimal check

    int n = (int)atoms.size();
    adj.assign(n, {});
    auto addBond = [&](const std::string& a, const std::string& b) {
        auto ia = nameIdx.find(a), ib = nameIdx.find(b);
        if (ia != nameIdx.end() && ib != nameIdx.end()) {
            adj[ia->second].push_back(ib->second);
            adj[ib->second].push_back(ia->second);
        }
    };
    for (const auto& b : bbBonds) addBond(b.first, b.second);
    for (const auto& b : scBonds) addBond(b.first, b.second);

    // Add bonds from H atoms to their parent heavy atoms.
    // H naming convention: H* is bonded to the heavy atom whose name
    // matches the part after 'H'.  E.g. HB2→CB, HG13→CG1, HA→CA, H→N.
    for (int i = 0; i < n; ++i) {
        const std::string& nm = atoms[i].name;
        if (nm.empty() || nm[0] != 'H') continue;
        if (nm == "X_PI") continue;
        // Determine parent heavy atom name
        std::string parent;
        if (nm == "H") {
            parent = "N";  // backbone amide H
        } else if (nm.size() >= 2) {
            std::string suffix = nm.substr(1); // everything after 'H'
            parent = "C" + suffix;             // try C+suffix first (most common)
            while (!parent.empty() && nameIdx.find(parent) == nameIdx.end()) {
                if (parent.size() > 1) parent.pop_back();
                else break;
            }
            if (nameIdx.find(parent) == nameIdx.end()) {
                // Try N, O, S prefixes for polar H
                for (const char* prefix : {"N", "O", "S"}) {
                    parent = std::string(prefix) + suffix;
                    while (!parent.empty() && nameIdx.find(parent) == nameIdx.end()) {
                        if (parent.size() > 1) parent.pop_back();
                        else break;
                    }
                    if (nameIdx.find(parent) != nameIdx.end()) break;
                }
            }
        }
        if (!parent.empty() && nameIdx.find(parent) != nameIdx.end()) {
            adj[i].push_back(nameIdx[parent]);
            adj[nameIdx[parent]].push_back(i);
        }
    }
    return true;
}

// ── Virtual hydrogen placement ──────────────────────────────────────────
// Add hydrogen atoms to a residue's FxAtom list for clash evaluation.
// FoldX uses H_vdw_radius = 0.99 Å for all H atoms.
// We place H using simple tetrahedral/trigonal geometry from heavy-atom
// positions (same principle as FoldX's LookupTransform, but computed
// on the fly from the current heavy-atom coordinates).
//
// H atoms are essential for:
//   - Inter-residue clash: H-X contacts dominate backbone clashes
//   - Intra-residue (torsional) clash: H-H contacts between adjacent
//     rotatable bonds are the primary signal for eclipsed conformations

/// Standard FxAtomParam for non-polar H (C-H bonds)
/// volume=0: H atoms don't contribute to Facc burial calculation.
/// vdw=0, enerG=0: H atoms don't contribute to VdW/solvation energy.
/// min_radius=0.99: FoldX H_vdw_radius for clash detection.
inline FxAtomParam HParamApolar() {
    return {0.0, 375, 375, 0.0, 0.0, 0.99, 0.99, 0, 0, 0};
}
/// Standard FxAtomParam for polar H (N-H, O-H, S-H bonds)
inline FxAtomParam HParamPolar() {
    return {0.0, 375, 375, 0.0, 0.0, 0.99, 0.99, 1, 0, 0};
}

/// Helper: normalize a vector, returns zero vector if degenerate
inline Vector3D Normalize(const Vector3D& v) {
    double len = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (len < 1e-8) return Vector3D(0,0,0);
    return Vector3D(v.x/len, v.y/len, v.z/len);
}

/// Helper: cross product
inline Vector3D Cross(const Vector3D& a, const Vector3D& b) {
    return Vector3D(a.y*b.z - a.z*b.y,
                    a.z*b.x - a.x*b.z,
                    a.x*b.y - a.y*b.x);
}

/// Place 1 H on a sp3 center with 3 heavy neighbors (CH1: Val CB, Ile CB, Thr CB, Leu CG)
inline Vector3D PlaceCH1(const Vector3D& center,
                         const Vector3D& a, const Vector3D& b, const Vector3D& c,
                         double bondLen = 1.09) {
    Vector3D da = Normalize(a - center);
    Vector3D db = Normalize(b - center);
    Vector3D dc = Normalize(c - center);
    Vector3D sum(da.x+db.x+dc.x, da.y+db.y+dc.y, da.z+db.z+dc.z);
    Vector3D dir = Normalize(sum);
    return Vector3D(center.x - bondLen*dir.x, center.y - bondLen*dir.y, center.z - bondLen*dir.z);
}

/// Place 2 H on a sp3 center with 2 heavy neighbors (CH2: most CB, CG, CD, CE)
inline void PlaceCH2(const Vector3D& center,
                     const Vector3D& a, const Vector3D& b,
                     Vector3D& h1, Vector3D& h2,
                     double bondLen = 1.09) {
    Vector3D da = Normalize(a - center);
    Vector3D db = Normalize(b - center);
    // Bisector of the two heavy neighbors
    Vector3D bisector(da.x+db.x, da.y+db.y, da.z+db.z);
    bisector = Normalize(bisector);
    // Normal to the a-center-b plane
    Vector3D normal = Normalize(Cross(da, db));
    // The two H lie in the plane that bisects the a-C-b angle and is
    // perpendicular to the a-C-b plane, on the side AWAY from the heavy
    // neighbours (the -bisector direction), symmetric about ±normal.  For
    // ideal tetrahedral geometry the H sit at ±54.735° from -bisector:
    //   cos(54.735°) = 0.57735 along -bisector,  sin = 0.81650 along ±normal.
    // (The earlier 0.3338/0.9426 along +bisector gave an H-C-H angle of 141°
    //  and placed each H ~1.0 Å off; see PlaceCH2 history.)
    double axial   = 0.57735 * bondLen;   // along -bisector (away from a,b)
    double lateral = 0.81650 * bondLen;   // perpendicular, in ±normal
    h1 = Vector3D(center.x - axial*bisector.x + lateral*normal.x,
                   center.y - axial*bisector.y + lateral*normal.y,
                   center.z - axial*bisector.z + lateral*normal.z);
    h2 = Vector3D(center.x - axial*bisector.x - lateral*normal.x,
                   center.y - axial*bisector.y - lateral*normal.y,
                   center.z - axial*bisector.z - lateral*normal.z);
}

/// Place 3 H on a sp3 center with 1 heavy neighbor (CH3: Ala CB, methyl groups)
inline void PlaceCH3(const Vector3D& center, const Vector3D& a,
                     const Vector3D& ref,   // any other atom for reference frame
                     Vector3D& h1, Vector3D& h2, Vector3D& h3,
                     double bondLen = 1.09) {
    Vector3D da = Normalize(a - center);
    // Build a perpendicular reference
    Vector3D dr = Normalize(ref - center);
    Vector3D perp1 = Normalize(Cross(da, dr));
    if (perp1.x*perp1.x + perp1.y*perp1.y + perp1.z*perp1.z < 0.01) {
        // da and dr are parallel, use arbitrary perpendicular
        Vector3D arb(1, 0, 0);
        if (std::fabs(da.x) > 0.9) arb = Vector3D(0, 1, 0);
        perp1 = Normalize(Cross(da, arb));
    }
    Vector3D perp2 = Cross(da, perp1);
    // Tetrahedral: angle from -da to each H is 109.5°
    // cos(109.5°) = -1/3, sin(109.5°) = sqrt(8/9) ≈ 0.9428
    double axial   = -0.3338 * bondLen;   // along -da
    double lateral =  0.9426 * bondLen;   // perpendicular
    // Three H atoms at 120° intervals around -da axis
    for (int k = 0; k < 3; ++k) {
        double angle = k * 2.0 * M_PI / 3.0;
        double px = lateral * (std::cos(angle) * perp1.x + std::sin(angle) * perp2.x);
        double py = lateral * (std::cos(angle) * perp1.y + std::sin(angle) * perp2.y);
        double pz = lateral * (std::cos(angle) * perp1.z + std::sin(angle) * perp2.z);
        Vector3D h(center.x + axial*(-da.x) + px,
                   center.y + axial*(-da.y) + py,
                   center.z + axial*(-da.z) + pz);
        if (k == 0) h1 = h;
        else if (k == 1) h2 = h;
        else h3 = h;
    }
}

/// Place 1 H on a trigonal sp2 center with 2 heavy neighbors (aromatic CH, amide NH)
inline Vector3D PlaceSP2H(const Vector3D& center,
                          const Vector3D& a, const Vector3D& b,
                          double bondLen = 1.08) {
    Vector3D da = Normalize(a - center);
    Vector3D db = Normalize(b - center);
    // H is opposite to the bisector of a and b, in the sp2 plane
    Vector3D bisector(da.x+db.x, da.y+db.y, da.z+db.z);
    Vector3D dir = Normalize(bisector);
    return Vector3D(center.x - bondLen*dir.x, center.y - bondLen*dir.y, center.z - bondLen*dir.z);
}

/// Place 2 H on a sp2 nitrogen with 1 heavy neighbor (NH2: Asn ND2, Gln NE2, Arg NH1/NH2)
inline void PlaceNH2(const Vector3D& center, const Vector3D& heavy,
                     const Vector3D& planeRef,  // atom to define the sp2 plane
                     Vector3D& h1, Vector3D& h2,
                     double bondLen = 1.01) {
    Vector3D dc = Normalize(heavy - center);
    Vector3D dr = Normalize(planeRef - center);
    // sp2 plane normal
    Vector3D normal = Normalize(Cross(dc, dr));
    // Two H at 120° from the heavy neighbor, in the sp2 plane
    // Rotate dc by ±120° around normal
    static constexpr double cos120 = -0.5, sin120 = 0.866025;
    // Perpendicular to dc in the sp2 plane
    Vector3D perp = Cross(normal, dc);
    h1 = Vector3D(center.x + bondLen*(cos120*dc.x + sin120*perp.x),
                  center.y + bondLen*(cos120*dc.y + sin120*perp.y),
                  center.z + bondLen*(cos120*dc.z + sin120*perp.z));
    h2 = Vector3D(center.x + bondLen*(cos120*dc.x - sin120*perp.x),
                  center.y + bondLen*(cos120*dc.y - sin120*perp.y),
                  center.z + bondLen*(cos120*dc.z - sin120*perp.z));
}

/// Add virtual hydrogen atoms to a residue's FxAtom list.
/// Call this after AnnotateHostFx, before scoring.
inline void AddVirtualHydrogens(std::vector<FxAtom>& atoms, char aa,
                                const std::string& res3) {
    // Build name→position map
    std::map<std::string, Vector3D> pos;
    for (const auto& a : atoms) {
        if (a.name == "X_PI") continue;
        pos[a.name] = a.pos;
    }

    auto has = [&](const char* n) { return pos.find(n) != pos.end(); };
    auto P   = [&](const char* n) -> Vector3D { return pos.at(n); };

    // Get resSeqIdx from the first real atom
    int resIdx = -1;
    for (const auto& a : atoms) {
        if (a.resSeqIdx >= 0) { resIdx = a.resSeqIdx; break; }
    }

    // Helper: append an H atom
    auto addH = [&](const std::string& name, const Vector3D& p, bool polar) {
        FxAtom h;
        h.name = name;
        h.pos = p;
        h.param = polar ? HParamPolar() : HParamApolar();
        h.charge = 0.0;
        h.resSeqIdx = resIdx;
        atoms.push_back(h);
    };

    // ── Backbone H ──
    // HA on CA (not GLY, not PRO)
    if (aa != 'G' && aa != 'P' && has("CA") && has("N") && has("C") && has("CB")) {
        Vector3D ha = PlaceCH1(P("CA"), P("N"), P("C"), P("CB"));
        addH("HA", ha, false);
    }
    // GLY: HA2, HA3 (two H on CA)
    if (aa == 'G' && has("CA") && has("N") && has("C")) {
        Vector3D h1, h2;
        PlaceCH2(P("CA"), P("N"), P("C"), h1, h2);
        addH("HA2", h1, false);
        addH("HA3", h2, false);
    }
    // HN on backbone N (not PRO)
    // HN is placed using CA and C as references, trans to C around N-CA
    if (aa != 'P' && has("N") && has("CA") && has("C")) {
        // N is sp2-like in peptide bond: planar with CA and previous C
        // Place HN opposite to the CA-N-C bisector
        Vector3D hn = PlaceSP2H(P("N"), P("CA"), P("C"), 1.01);
        addH("H", hn, true);
    }

    // ── Sidechain H (amino-acid-specific) ──

    // Helper for common patterns
    auto addCH2 = [&](const char* center, const char* nbr1, const char* nbr2,
                       const std::string& hname1, const std::string& hname2) {
        if (has(center) && has(nbr1) && has(nbr2)) {
            Vector3D h1, h2;
            PlaceCH2(P(center), P(nbr1), P(nbr2), h1, h2);
            addH(hname1, h1, false);
            addH(hname2, h2, false);
        }
    };

    auto addCH3 = [&](const char* center, const char* heavy, const char* ref,
                       const std::string& h1n, const std::string& h2n, const std::string& h3n) {
        if (has(center) && has(heavy) && has(ref)) {
            Vector3D h1, h2, h3;
            PlaceCH3(P(center), P(heavy), P(ref), h1, h2, h3);
            addH(h1n, h1, false);
            addH(h2n, h2, false);
            addH(h3n, h3, false);
        }
    };

    switch (aa) {
    case 'A':
        // CB is CH3 (bonded to CA only)
        addCH3("CB", "CA", "N", "HB1", "HB2", "HB3");
        break;

    case 'V':
        // CB is CH1 (CA, CG1, CG2)
        if (has("CB") && has("CA") && has("CG1") && has("CG2")) {
            addH("HB", PlaceCH1(P("CB"), P("CA"), P("CG1"), P("CG2")), false);
        }
        addCH3("CG1", "CB", "CA", "HG11", "HG12", "HG13");
        addCH3("CG2", "CB", "CA", "HG21", "HG22", "HG23");
        break;

    case 'I':
        // CB is CH1 (CA, CG1, CG2)
        if (has("CB") && has("CA") && has("CG1") && has("CG2")) {
            addH("HB", PlaceCH1(P("CB"), P("CA"), P("CG1"), P("CG2")), false);
        }
        addCH2("CG1", "CB", "CD1", "HG12", "HG13");
        addCH3("CG2", "CB", "CA", "HG21", "HG22", "HG23");
        addCH3("CD1", "CG1", "CB", "HD11", "HD12", "HD13");
        break;

    case 'L':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        // CG is CH1 (CB, CD1, CD2)
        if (has("CG") && has("CB") && has("CD1") && has("CD2")) {
            addH("HG", PlaceCH1(P("CG"), P("CB"), P("CD1"), P("CD2")), false);
        }
        addCH3("CD1", "CG", "CB", "HD11", "HD12", "HD13");
        addCH3("CD2", "CG", "CB", "HD21", "HD22", "HD23");
        break;

    case 'S':
        addCH2("CB", "CA", "OG", "HB2", "HB3");
        // HG on OG (hydroxyl)
        if (has("OG") && has("CB") && has("CA")) {
            addH("HG", PlaceSP2H(P("OG"), P("CB"), P("CA"), 0.96), true);
        }
        break;

    case 'T':
        // CB is CH1 (CA, OG1, CG2)
        if (has("CB") && has("CA") && has("OG1") && has("CG2")) {
            addH("HB", PlaceCH1(P("CB"), P("CA"), P("OG1"), P("CG2")), false);
        }
        // HG1 on OG1
        if (has("OG1") && has("CB") && has("CA")) {
            addH("HG1", PlaceSP2H(P("OG1"), P("CB"), P("CA"), 0.96), true);
        }
        addCH3("CG2", "CB", "CA", "HG21", "HG22", "HG23");
        break;

    case 'C':
        addCH2("CB", "CA", "SG", "HB2", "HB3");
        // HG on SG (thiol)
        if (has("SG") && has("CB") && has("CA")) {
            addH("HG", PlaceSP2H(P("SG"), P("CB"), P("CA"), 1.34), true);
        }
        break;

    case 'M':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        addCH2("CG", "CB", "SD", "HG2", "HG3");
        addCH3("CE", "SD", "CG", "HE1", "HE2", "HE3");
        break;

    case 'F':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        // Aromatic H: sp2, each ring C with 2 ring neighbors gets 1 H
        if (has("CD1") && has("CG") && has("CE1"))
            addH("HD1", PlaceSP2H(P("CD1"), P("CG"), P("CE1"), 1.08), false);
        if (has("CD2") && has("CG") && has("CE2"))
            addH("HD2", PlaceSP2H(P("CD2"), P("CG"), P("CE2"), 1.08), false);
        if (has("CE1") && has("CD1") && has("CZ"))
            addH("HE1", PlaceSP2H(P("CE1"), P("CD1"), P("CZ"), 1.08), false);
        if (has("CE2") && has("CD2") && has("CZ"))
            addH("HE2", PlaceSP2H(P("CE2"), P("CD2"), P("CZ"), 1.08), false);
        if (has("CZ") && has("CE1") && has("CE2"))
            addH("HZ", PlaceSP2H(P("CZ"), P("CE1"), P("CE2"), 1.08), false);
        break;

    case 'Y':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        if (has("CD1") && has("CG") && has("CE1"))
            addH("HD1", PlaceSP2H(P("CD1"), P("CG"), P("CE1"), 1.08), false);
        if (has("CD2") && has("CG") && has("CE2"))
            addH("HD2", PlaceSP2H(P("CD2"), P("CG"), P("CE2"), 1.08), false);
        if (has("CE1") && has("CD1") && has("CZ"))
            addH("HE1", PlaceSP2H(P("CE1"), P("CD1"), P("CZ"), 1.08), false);
        if (has("CE2") && has("CD2") && has("CZ"))
            addH("HE2", PlaceSP2H(P("CE2"), P("CD2"), P("CZ"), 1.08), false);
        // OH gets HH
        if (has("OH") && has("CZ") && has("CE1"))
            addH("HH", PlaceSP2H(P("OH"), P("CZ"), P("CE1"), 0.96), true);
        break;

    case 'W':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        // Indole ring H atoms
        if (has("CD1") && has("CG") && has("NE1"))
            addH("HD1", PlaceSP2H(P("CD1"), P("CG"), P("NE1"), 1.08), false);
        if (has("NE1") && has("CD1") && has("CE2"))
            addH("HE1", PlaceSP2H(P("NE1"), P("CD1"), P("CE2"), 1.01), true);
        if (has("CE3") && has("CD2") && has("CZ3"))
            addH("HE3", PlaceSP2H(P("CE3"), P("CD2"), P("CZ3"), 1.08), false);
        if (has("CZ2") && has("CE2") && has("CH2"))
            addH("HZ2", PlaceSP2H(P("CZ2"), P("CE2"), P("CH2"), 1.08), false);
        if (has("CZ3") && has("CE3") && has("CH2"))
            addH("HZ3", PlaceSP2H(P("CZ3"), P("CE3"), P("CH2"), 1.08), false);
        if (has("CH2") && has("CZ2") && has("CZ3"))
            addH("HH2", PlaceSP2H(P("CH2"), P("CZ2"), P("CZ3"), 1.08), false);
        break;

    case 'H':  // HIS: both ND1 and NE2 protonated (standard treatment)
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        if (has("ND1") && has("CG") && has("CE1"))
            addH("HD1", PlaceSP2H(P("ND1"), P("CG"), P("CE1"), 1.01), true);   // N-H donor
        if (has("CE1") && has("ND1") && has("NE2"))
            addH("HE1", PlaceSP2H(P("CE1"), P("ND1"), P("NE2"), 1.08), false); // ring C-H
        if (has("NE2") && has("CD2") && has("CE1"))
            addH("HE2", PlaceSP2H(P("NE2"), P("CD2"), P("CE1"), 1.01), true);  // N-H donor
        if (has("CD2") && has("CG") && has("NE2"))
            addH("HD2", PlaceSP2H(P("CD2"), P("CG"), P("NE2"), 1.08), false);  // ring C-H
        break;

    case 'D':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        // No H on CG (sp2 carboxylate), no H on OD1/OD2 (ionized)
        break;

    case 'N':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        // ND2 gets 2 H (amide NH2)
        if (has("ND2") && has("CG") && has("OD1")) {
            Vector3D h1, h2;
            PlaceNH2(P("ND2"), P("CG"), P("OD1"), h1, h2, 1.01);
            addH("HD21", h1, true);
            addH("HD22", h2, true);
        }
        break;

    case 'E':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        addCH2("CG", "CB", "CD", "HG2", "HG3");
        // No H on CD (sp2 carboxylate)
        break;

    case 'Q':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        addCH2("CG", "CB", "CD", "HG2", "HG3");
        // NE2 gets 2 H (amide NH2)
        if (has("NE2") && has("CD") && has("OE1")) {
            Vector3D h1, h2;
            PlaceNH2(P("NE2"), P("CD"), P("OE1"), h1, h2, 1.01);
            addH("HE21", h1, true);
            addH("HE22", h2, true);
        }
        break;

    case 'K':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        addCH2("CG", "CB", "CD", "HG2", "HG3");
        addCH2("CD", "CG", "CE", "HD2", "HD3");
        addCH2("CE", "CD", "NZ", "HE2", "HE3");
        // NZ: NH3+ (3 H), place as CH3-like
        if (has("NZ") && has("CE") && has("CD")) {
            Vector3D h1, h2, h3;
            PlaceCH3(P("NZ"), P("CE"), P("CD"), h1, h2, h3, 1.01);
            addH("HZ1", h1, true);
            addH("HZ2", h2, true);
            addH("HZ3", h3, true);
        }
        break;

    case 'R':
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        addCH2("CG", "CB", "CD", "HG2", "HG3");
        addCH2("CD", "CG", "NE", "HD2", "HD3");
        // NE gets 1 H (sp2)
        if (has("NE") && has("CD") && has("CZ"))
            addH("HE", PlaceSP2H(P("NE"), P("CD"), P("CZ"), 1.01), true);
        // NH1, NH2 each get 2 H (sp2 NH2)
        if (has("NH1") && has("CZ") && has("NE")) {
            Vector3D h1, h2;
            PlaceNH2(P("NH1"), P("CZ"), P("NE"), h1, h2, 1.01);
            addH("HH11", h1, true);
            addH("HH12", h2, true);
        }
        if (has("NH2") && has("CZ") && has("NH1")) {
            Vector3D h1, h2;
            PlaceNH2(P("NH2"), P("CZ"), P("NH1"), h1, h2, 1.01);
            addH("HH21", h1, true);
            addH("HH22", h2, true);
        }
        break;

    case 'P':
        // PRO: CB-CG-CD are all CH2; no HN
        addCH2("CB", "CA", "CG", "HB2", "HB3");
        addCH2("CG", "CB", "CD", "HG2", "HG3");
        addCH2("CD", "CG", "N",  "HD2", "HD3");
        break;

    default:
        break;
    }
}

// ── Atom parameter lookup ───────────────────────────────────────────
// Key: "RES:ATOM" e.g. "ALA:CB", "PHE:CE1". Per-atom volume, Occ and Occmax
// of the FoldX 2014 parameter set, read from the compiled parameter library
// (atomparams.h, libatomparams.a). Only these three fields are set; the other
// FxAtomParam fields are zero. Unknown atoms get one shared fallback object.

inline const FxAtomParam& FxLookup(const std::string& res3, const std::string& atomName)
{
    // Static lookup map built on first call
    static const auto& table = *[]() {
        auto* m = new std::map<std::string, FxAtomParam>();
        const int n = ap_n_atoms();
        for (int i = 0; i < n; ++i) {
            const char* res = nullptr; const char* atm = nullptr;
            FxAtomParam p{};
            if (ap_atom(i, &res, &atm, &p.volume, &p.Occ, &p.Occmax))
                (*m)[std::string(res) + ":" + std::string(atm)] = p;
        }
        return m;
    }();

    std::string key = res3 + ":" + atomName;
    auto it = table.find(key);
    if (it != table.end()) return it->second;

    static const FxAtomParam fallback = [] {
        FxAtomParam p{};
        ap_fallback(&p.volume, &p.Occ, &p.Occmax);
        return p;
    }();
    return fallback;
}


/// Occupancy (FoldX "contrat") of a probe position: the summed volume of the
/// environment atoms within cutoff, with FoldX's Gaussian distance decay
/// volume * exp(-d^2 / 24.5). Atoms closer than 0.1 A (the probe itself) are
/// skipped. Split out of ComputeFacc so callers that need the raw occupancy
/// (e.g. FoldX-style per-residue burial as a ratio of sums) can get it.
inline double ComputeOccupancy(const Vector3D& probePos,
                               const std::vector<FxAtom>& envAtoms,
                               double cutoff = 6.0)
{
    double cut2 = cutoff * cutoff;
    double contrat = 0.0;

    for (const auto& e : envAtoms) {
        double d2 = Dist2(probePos, e.pos);
        if (d2 > cut2 || d2 < 0.01) continue;
        // FoldX uses Gaussian distance decay: volume * exp(-d^2 / 24.5)
        contrat += e.param.volume * std::exp(-d2 / 24.5);
    }
    return contrat;
}

} // namespace ZipperScore
