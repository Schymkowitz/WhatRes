// BondGeometry.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Extracts covalent bond lengths and backbone bond angles from PDB atom
// coordinates.  These are NOT available from FoldX SequenceDetail, so this
// fills the gap for the structure survey.
//
// Per-residue output:
//   bond_N_CA      N–Cα bond length (Å)     expected ~1.47
//   bond_CA_C      Cα–C bond length (Å)      expected ~1.52
//   bond_C_N       C–N peptide bond to next residue (Å)  expected ~1.33
//   bond_C_O       C=O bond length (Å)       expected ~1.24
//   bond_CA_CB     Cα–Cβ bond length (Å)     expected ~1.52 (absent for Gly)
//   angle_N_CA_C   N–Cα–C bond angle (°)     expected ~111
//   angle_CA_C_N   Cα–C–N bond angle (°)     expected ~117
//   angle_C_N_CA   C–N–Cα bond angle (°)     expected ~121
//
// Usage:
//   vector<BondGeometryResult> results = BondGeometry::Compute(peptide);
//
// Designed for speed: pure coordinate math, no allocations beyond output.
// -----------------------------------------------------------------------------

#pragma once
#include <vector>
#include <string>
#include <cmath>
#include "Peptide.h"
#include "Residue.h"
#include "Atom.h"

struct BondGeometryResult
{
    int    resIndex     = -1;     // index in Peptide::reslist
    char   chain        = ' ';
    char   aa1          = '?';
    int    resnum       = 0;

    // Bond lengths (Å).  NaN if atoms not found.
    double bond_N_CA    = NAN;
    double bond_CA_C    = NAN;
    double bond_C_N     = NAN;    // peptide bond to NEXT residue's N
    double bond_C_O     = NAN;
    double bond_CA_CB   = NAN;    // NaN for Gly

    // Bond angles (degrees).  NaN if atoms not found.
    double angle_N_CA_C = NAN;
    double angle_CA_C_N = NAN;    // to NEXT residue's N
    double angle_C_N_CA = NAN;    // from PREVIOUS residue's C

    // Cβ deviation (Å): distance of the observed Cβ from its ideal position
    // reconstructed from N, Cα, C (Lovell et al. 2003).  NaN for Gly / missing.
    double cb_deviation = NAN;

    // Sidechain dihedrals (degrees).  NaN where the residue type has no such
    // dihedral or atoms are missing.  chi1 = N-CA-CB-Xγ, chi2 = CA-CB-Xγ-Xδ,
    // chi3 = CB-Xγ-Xδ-Xε, chi4 = Xγ-Xδ-Xε-Xζ.
    double chi1 = NAN;
    double chi2 = NAN;
    double chi3 = NAN;
    double chi4 = NAN;

    // Peptide-bond planarity: |ω| deviation from the nearest planar value
    // (0 = cis or trans, larger = twisted), degrees [0,90].  Uses next residue.
    double peptide_planarity = NAN;

    // Rotamer outlier proxy: angular distance of χ1 from the nearest staggered
    // well (-60/+60/180), degrees [0,~60].  NaN for Gly/Ala/Pro / missing χ1.
    double rotamer_outlier = NAN;

    // CaBLAM Cα-CO virtual-dihedral measures (Williams et al. 2018), backbone /
    // low-resolution validation, robust where refinement restrains other metrics.
    //   cablam_mu_in  = dihedral CA(i-2)-CA(i-1)-CA(i)-CA(i+1)   (deg)
    //   cablam_mu_out = dihedral CA(i-1)-CA(i)-CA(i+1)-CA(i+2)   (deg)
    //   cablam_nu     = CO peptide-plane dihedral from O(i-1),O(i) vs CA trace (deg)
    //   cablam_ca_virtual = virtual angle CA(i-1)-CA(i)-CA(i+1)  (deg)
    double cablam_mu_in      = NAN;
    double cablam_mu_out     = NAN;
    double cablam_nu         = NAN;
    double cablam_ca_virtual = NAN;

    // cis-peptide indicator: 1.0 if the C(i)–N(i+1) peptide bond is cis (|ω|<30°),
    // else 0.0.  NaN if ω unavailable (e.g. chain end).
    double cis_peptide = NAN;

    // Cα chirality: improper dihedral N-CA-C-CB (deg).  L-amino acids cluster
    // tightly (empirically ≈ -122° for this atom ordering); a D-residue or
    // chirality violation shifts it by ~120°.  NaN for Gly / missing atoms.
    double ca_chirality = NAN;
};

class BondGeometry
{
public:
    // Compute bond geometry for all residues in a chain.
    static std::vector<BondGeometryResult> Compute(const Peptide& pep)
    {
        std::vector<BondGeometryResult> results;
        int n = (int)pep.reslist.size();
        if (n == 0) return results;
        results.resize(n);

        for (int i = 0; i < n; i++)
        {
            const Residue& r = pep.reslist[i];
            BondGeometryResult& bg = results[i];
            bg.resIndex = i;
            bg.chain    = r.chain;
            bg.aa1      = r.OneLetter();
            bg.resnum   = r.number_in_pdb;

            // Get backbone atom positions
            double Nx, Ny, Nz;     bool hasN  = findAtom(r, "N",  Nx, Ny, Nz);
            double CAx, CAy, CAz;  bool hasCA = findAtom(r, "CA", CAx, CAy, CAz);
            double Cx, Cy, Cz;     bool hasC  = findAtom(r, "C",  Cx, Cy, Cz);
            double Ox, Oy, Oz;     bool hasO  = findAtom(r, "O",  Ox, Oy, Oz);
            double CBx, CBy, CBz;  bool hasCB = findAtom(r, "CB", CBx, CBy, CBz);

            // Intra-residue bond lengths
            if (hasN  && hasCA) bg.bond_N_CA  = dist(Nx,Ny,Nz, CAx,CAy,CAz);
            if (hasCA && hasC)  bg.bond_CA_C  = dist(CAx,CAy,CAz, Cx,Cy,Cz);
            if (hasC  && hasO)  bg.bond_C_O   = dist(Cx,Cy,Cz, Ox,Oy,Oz);
            if (hasCA && hasCB) bg.bond_CA_CB = dist(CAx,CAy,CAz, CBx,CBy,CBz);

            // Cβ deviation from ideal position (needs full backbone + Cβ)
            if (hasN && hasCA && hasC && hasCB)
                bg.cb_deviation = cbDeviation(Nx,Ny,Nz, CAx,CAy,CAz,
                                              Cx,Cy,Cz, CBx,CBy,CBz);

            // Sidechain dihedrals chi1..chi4 (per amino-acid type)
            computeChi(r, bg.aa1, bg.chi1, bg.chi2, bg.chi3, bg.chi4);
            bg.rotamer_outlier = rotamerOutlier(bg.aa1, bg.chi1);

            // Cα chirality improper N-CA-C-CB (L ≈ +34°)
            if (hasN && hasCA && hasC && hasCB)
                bg.ca_chirality = dihedral(Nx,Ny,Nz, CAx,CAy,CAz,
                                           Cx,Cy,Cz, CBx,CBy,CBz);

            // CaBLAM Cα/CO virtual-dihedral measures (needs i-2..i+2 neighbours)
            cablam(pep, i, bg.cablam_mu_in, bg.cablam_mu_out,
                   bg.cablam_nu, bg.cablam_ca_virtual);

            // Intra-residue angle: N-CA-C
            if (hasN && hasCA && hasC)
                bg.angle_N_CA_C = angle(Nx,Ny,Nz, CAx,CAy,CAz, Cx,Cy,Cz);

            // Inter-residue: peptide bond C(i)–N(i+1) and angles
            if (i + 1 < n && hasC)
            {
                const Residue& rNext = pep.reslist[i + 1];
                double Nx2, Ny2, Nz2;
                double CAx2, CAy2, CAz2;
                bool hasN2  = findAtom(rNext, "N",  Nx2, Ny2, Nz2);
                bool hasCA2 = findAtom(rNext, "CA", CAx2, CAy2, CAz2);

                if (hasN2)
                    bg.bond_C_N = dist(Cx,Cy,Cz, Nx2,Ny2,Nz2);

                if (hasCA && hasN2)
                    bg.angle_CA_C_N = angle(CAx,CAy,CAz, Cx,Cy,Cz, Nx2,Ny2,Nz2);

                if (hasN2 && hasCA2)
                    results[i + 1].angle_C_N_CA = angle(Cx,Cy,Cz, Nx2,Ny2,Nz2, CAx2,CAy2,CAz2);

                // Peptide-bond planarity from ω = Cα(i)-C(i)-N(i+1)-Cα(i+1)
                if (hasCA && hasN2 && hasCA2) {
                    double omega = dihedral(CAx,CAy,CAz, Cx,Cy,Cz,
                                            Nx2,Ny2,Nz2, CAx2,CAy2,CAz2);
                    if (!std::isnan(omega)) {
                        double ao = std::fabs(omega);            // 0=cis, 180=trans
                        bg.peptide_planarity = std::min(std::fabs(ao - 180.0), ao);
                        bg.cis_peptide = (ao < 30.0) ? 1.0 : 0.0;
                    }
                }
            }
        }

        return results;
    }

    // Header string for TSV output
    static const char* Header()
    {
        return "bond_N_CA\tbond_CA_C\tbond_C_N\tbond_C_O\tbond_CA_CB"
               "\tangle_N_CA_C\tangle_CA_C_N\tangle_C_N_CA\tcb_deviation"
               "\tchi1\tchi2\tchi3\tchi4\tpeptide_planarity\trotamer_outlier"
               "\tcablam_mu_in\tcablam_mu_out\tcablam_nu\tcablam_ca_virtual"
               "\tcis_peptide\tca_chirality";
    }

    // Format one result as TSV columns (matches Header)
    static std::string ToTSV(const BondGeometryResult& bg)
    {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f\t%.1f\t%.1f\t%.3f"
                 "\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f"
                 "\t%.1f\t%.1f\t%.1f\t%.1f\t%.0f\t%.1f",
                 bg.bond_N_CA, bg.bond_CA_C, bg.bond_C_N,
                 bg.bond_C_O, bg.bond_CA_CB,
                 bg.angle_N_CA_C, bg.angle_CA_C_N, bg.angle_C_N_CA,
                 bg.cb_deviation, bg.chi1, bg.chi2, bg.chi3, bg.chi4,
                 bg.peptide_planarity, bg.rotamer_outlier,
                 bg.cablam_mu_in, bg.cablam_mu_out, bg.cablam_nu,
                 bg.cablam_ca_virtual, bg.cis_peptide, bg.ca_chirality);
        // Replace "nan" with empty for cleaner output
        std::string s(buf);
        for (size_t pos = 0; (pos = s.find("nan", pos)) != std::string::npos; )
            s.replace(pos, 3, "");
        return s;
    }

private:
    static bool findAtom(const Residue& r, const char* name,
                         double& ox, double& oy, double& oz)
    {
        for (const auto& a : r.atoms)
        {
            // Exact match on the trimmed PDB atom name.  Substring matching is
            // wrong: querying "C" also matches "CA"/"CB"/"CG"..., which silently
            // corrupted bond_CA_C (always 0), bond_C_N, bond_C_O and the
            // backbone angles.  Atom names are stored trimmed, but trim defensively.
            const std::string& an = a.atom_name;
            size_t b = an.find_first_not_of(' ');
            if (b != std::string::npos)
            {
                size_t e = an.find_last_not_of(' ');
                if (an.compare(b, e - b + 1, name) == 0)
                {
                    ox = a.x; oy = a.y; oz = a.z;
                    return true;
                }
            }
        }
        return false;
    }

    static double dist(double x1, double y1, double z1,
                       double x2, double y2, double z2)
    {
        double dx = x1-x2, dy = y1-y2, dz = z1-z2;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    static double angle(double x1, double y1, double z1,   // atom A
                        double x2, double y2, double z2,    // atom B (vertex)
                        double x3, double y3, double z3)    // atom C
    {
        // Angle A-B-C at vertex B
        double bax = x1-x2, bay = y1-y2, baz = z1-z2;
        double bcx = x3-x2, bcy = y3-y2, bcz = z3-z2;
        double dot  = bax*bcx + bay*bcy + baz*bcz;
        double magA = std::sqrt(bax*bax + bay*bay + baz*baz);
        double magC = std::sqrt(bcx*bcx + bcy*bcy + bcz*bcz);
        if (magA < 1e-9 || magC < 1e-9) return NAN;
        double cosTheta = dot / (magA * magC);
        // Clamp to [-1,1] for numerical safety
        if (cosTheta >  1.0) cosTheta =  1.0;
        if (cosTheta < -1.0) cosTheta = -1.0;
        return std::acos(cosTheta) * 180.0 / 3.14159265358979323846;
    }

    // Distance of the observed Cβ from its ideal position reconstructed from
    // N, Cα, C.  Geometry from Lovell et al. (2003), "Structure validation by
    // Cα geometry": Cα–Cβ = 1.55 Å, 51.1° out of the N–Cα–C plane, 1.81° from
    // the C–Cα–N angle bisector.  L-amino-acid chirality: normal = (N−Cα)×(C−Cα).
    static double cbDeviation(double Nx,double Ny,double Nz,
                              double CAx,double CAy,double CAz,
                              double Cx,double Cy,double Cz,
                              double CBx,double CBy,double CBz)
    {
        const double PI = 3.14159265358979323846;
        auto unit = [](double& x,double& y,double& z) {
            double L = std::sqrt(x*x + y*y + z*z);
            if (L > 1e-9) { x/=L; y/=L; z/=L; }
        };
        double nx=Nx-CAx, ny=Ny-CAy, nz=Nz-CAz; unit(nx,ny,nz);   // CA→N
        double cx=Cx-CAx, cy=Cy-CAy, cz=Cz-CAz; unit(cx,cy,cz);   // CA→C
        double bx=nx+cx, by=ny+cy, bz=nz+cz;     unit(bx,by,bz);   // bisector
        // normal to the N–CA–C plane (L-chirality)
        double mx=ny*cz-nz*cy, my=nz*cx-nx*cz, mz=nx*cy-ny*cx; unit(mx,my,mz);
        // in-plane axis perpendicular to the bisector
        double px=my*bz-mz*by, py=mz*bx-mx*bz, pz=mx*by-my*bx;
        const double tb = 1.81 * PI/180.0, tp = 51.1 * PI/180.0;
        double ix = -std::cos(tb)*bx + std::sin(tb)*px;
        double iy = -std::cos(tb)*by + std::sin(tb)*py;
        double iz = -std::cos(tb)*bz + std::sin(tb)*pz;
        double dx = std::cos(tp)*ix + std::sin(tp)*mx;
        double dy = std::cos(tp)*iy + std::sin(tp)*my;
        double dz = std::cos(tp)*iz + std::sin(tp)*mz;
        double cbx = CAx + 1.55*dx, cby = CAy + 1.55*dy, cbz = CAz + 1.55*dz;
        return dist(cbx,cby,cbz, CBx,CBy,CBz);
    }

    // Dihedral A-B-C-D in degrees, range (-180, 180].
    static double dihedral(double ax,double ay,double az, double bx,double by,double bz,
                           double cx,double cy,double cz, double dx,double dy,double dz)
    {
        double b1x=bx-ax, b1y=by-ay, b1z=bz-az;
        double b2x=cx-bx, b2y=cy-by, b2z=cz-bz;
        double b3x=dx-cx, b3y=dy-cy, b3z=dz-cz;
        // n1 = b1 × b2, n2 = b2 × b3
        double n1x=b1y*b2z-b1z*b2y, n1y=b1z*b2x-b1x*b2z, n1z=b1x*b2y-b1y*b2x;
        double n2x=b2y*b3z-b2z*b3y, n2y=b2z*b3x-b2x*b3z, n2z=b2x*b3y-b2y*b3x;
        double l1=std::sqrt(n1x*n1x+n1y*n1y+n1z*n1z);
        double l2=std::sqrt(n2x*n2x+n2y*n2y+n2z*n2z);
        if (l1<1e-9 || l2<1e-9) return NAN;
        double cosA=(n1x*n2x+n1y*n2y+n1z*n2z)/(l1*l2);
        if (cosA> 1.0) cosA= 1.0;
        if (cosA<-1.0) cosA=-1.0;
        double ang=std::acos(cosA)*180.0/3.14159265358979323846;
        // sign from (n1 × n2) · b2
        double cxv=n1y*n2z-n1z*n2y, cyv=n1z*n2x-n1x*n2z, czv=n1x*n2y-n1y*n2x;
        if (cxv*b2x+cyv*b2y+czv*b2z < 0.0) ang=-ang;
        return ang;
    }

    // Angular distance of χ1 from the nearest staggered rotamer well
    // (-60 / +60 / 180), in degrees [0,~60].  NaN for residues without a
    // meaningful staggered χ1 (Gly, Ala, ring-constrained Pro) or missing χ1.
    static double rotamerOutlier(char aa1, double chi1)
    {
        if (std::isnan(chi1) || aa1 == 'G' || aa1 == 'A' || aa1 == 'P') return NAN;
        const double wells[3] = {-60.0, 60.0, 180.0};
        double best = 1e9;
        for (double w : wells) {
            double d = std::fabs(chi1 - w);
            if (d > 180.0) d = 360.0 - d;   // wrap
            if (d < best) best = d;
        }
        return best;
    }

    // Perpendicular foot of point b1 onto the line a1-a2 (cctbx cablam perptersect).
    static void perptersect(const double a1[3], const double a2[3],
                            const double b1[3], double out[3])
    {
        double A[3] = {a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]};
        double denom = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
        if (denom < 1e-12) denom = 1e-12;
        double t = (A[0]*(b1[0]-a1[0]) + A[1]*(b1[1]-a1[1]) + A[2]*(b1[2]-a1[2])) / denom;
        out[0] = a1[0] + A[0]*t; out[1] = a1[1] + A[1]*t; out[2] = a1[2] + A[2]*t;
    }

    // ---- CaBLAM: Ca-Based Low-resolution Annotation Method (MolProbity lineage) ----
    // Backbone-conformation measures that use ONLY the Ca trace and the carbonyl
    // oxygens, deliberately ignoring the exact peptide N/C positions. Those positions
    // are poorly determined at low resolution and in cryo-EM, whereas the Ca trace and
    // the carbonyl O stay reliable, so these measures flag implausible backbone
    // conformation in a way that is robust to resolution.
    //
    // Computed per residue i (needs Ca of i-2..i+2 and carbonyl O of i-1, i):
    //   mu_in   = dihedral CA(i-2)-CA(i-1)-CA(i)-CA(i+1)   Ca virtual torsion, incoming
    //   mu_out  = dihedral CA(i-1)-CA(i)-CA(i+1)-CA(i+2)   Ca virtual torsion, outgoing
    //   ca_virt = angle    CA(i-1)-CA(i)-CA(i+1)           Ca virtual bond angle
    //   nu      = carbonyl-orientation dihedral (see below)
    //
    // nu construction: project the carbonyl O of i-1 perpendicularly onto the line
    // CA(i-1)-CA(i), and the carbonyl O of i onto CA(i)-CA(i+1) (perptersect, the
    // perpendicular foot), giving X1 and X2; then nu = dihedral(O(i-1), X1, X2, O(i)).
    // This captures how successive carbonyls twist relative to the Ca trace without
    // trusting the precise O positions.
    //
    // We compute the raw geometric values only. We do NOT apply MolProbity's CaBLAM
    // contour lookup that turns them into an outlier flag; instead these values are
    // scored against our own high-resolution reference distributions downstream.
    // Math matches the cctbx cablam source (calculate_mu / calculate_nu /
    // calculate_ca_virtual_angle), with res0..res4 = i-2..i+2. Reference: MolProbity
    // CaBLAM (Williams et al. 2018; detailed in Prisant et al. 2020); confirm the
    // canonical citation before publishing.
    static void cablam(const Peptide& pep, int i,
                       double& mu_in, double& mu_out, double& nu, double& ca_virt)
    {
        mu_in = mu_out = nu = ca_virt = NAN;
        int n = (int)pep.reslist.size();
        double ca[5][3]; bool hca[5];
        for (int k = 0; k < 5; ++k) {
            int idx = i - 2 + k;
            hca[k] = (idx >= 0 && idx < n) &&
                     findAtom(pep.reslist[idx], "CA", ca[k][0], ca[k][1], ca[k][2]);
        }
        if (hca[0] && hca[1] && hca[2] && hca[3])
            mu_in  = dihedral(ca[0][0],ca[0][1],ca[0][2], ca[1][0],ca[1][1],ca[1][2],
                              ca[2][0],ca[2][1],ca[2][2], ca[3][0],ca[3][1],ca[3][2]);
        if (hca[1] && hca[2] && hca[3] && hca[4])
            mu_out = dihedral(ca[1][0],ca[1][1],ca[1][2], ca[2][0],ca[2][1],ca[2][2],
                              ca[3][0],ca[3][1],ca[3][2], ca[4][0],ca[4][1],ca[4][2]);
        if (hca[1] && hca[2] && hca[3])
            ca_virt = angle(ca[1][0],ca[1][1],ca[1][2], ca[2][0],ca[2][1],ca[2][2],
                            ca[3][0],ca[3][1],ca[3][2]);
        // nu: carbonyl O of i-1 and i, perpendicular-projected onto the CA trace
        double o1[3], o2[3];
        bool ho1 = (i-1 >= 0)     && findAtom(pep.reslist[i-1], "O", o1[0],o1[1],o1[2]);
        bool ho2 = (i >= 0 && i<n) && findAtom(pep.reslist[i],   "O", o2[0],o2[1],o2[2]);
        if (hca[1] && hca[2] && hca[3] && ho1 && ho2) {
            double X1[3], X2[3];
            perptersect(ca[1], ca[2], o1, X1);
            perptersect(ca[2], ca[3], o2, X2);
            nu = dihedral(o1[0],o1[1],o1[2], X1[0],X1[1],X1[2],
                          X2[0],X2[1],X2[2], o2[0],o2[1],o2[2]);
        }
    }

    // Dihedral from four atoms given by name; NaN if any is missing.
    static double chiByNames(const Residue& r, const char* a, const char* b,
                             const char* c, const char* d)
    {
        double ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz;
        if (!findAtom(r,a,ax,ay,az) || !findAtom(r,b,bx,by,bz) ||
            !findAtom(r,c,cx,cy,cz) || !findAtom(r,d,dx,dy,dz)) return NAN;
        return dihedral(ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz);
    }

    // Sidechain dihedrals chi1..chi4 (per amino-acid type).  Atom definitions are
    // the IUPAC/standard ones (chi1 = N-CA-CB-Xγ, etc.); chi1/chi2 mirror
    // AmyloidRotamerExtractor::Chi1GammaAtom / Chi2Atoms.  Pro chi2 is omitted to
    // match that convention.  NaN where a dihedral does not apply.
    static void computeChi(const Residue& r, char aa1,
                           double& c1, double& c2, double& c3, double& c4)
    {
        c1 = c2 = c3 = c4 = NAN;
        switch (aa1) {
            case 'S': c1 = chiByNames(r,"N","CA","CB","OG");  break;
            case 'T': c1 = chiByNames(r,"N","CA","CB","OG1"); break;
            case 'C': c1 = chiByNames(r,"N","CA","CB","SG");  break;
            case 'V': c1 = chiByNames(r,"N","CA","CB","CG1"); break;
            case 'P': c1 = chiByNames(r,"N","CA","CB","CG");  break;  // ring; no chi2
            case 'I':
                c1 = chiByNames(r,"N","CA","CB","CG1");
                c2 = chiByNames(r,"CA","CB","CG1","CD1"); break;
            case 'L':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","CD1"); break;
            case 'F': case 'Y': case 'W': case 'H':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","CD1"); break;
            case 'D':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","OD1"); break;
            case 'N':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","OD1"); break;
            case 'M':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","SD");
                c3 = chiByNames(r,"CB","CG","SD","CE"); break;
            case 'E':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","CD");
                c3 = chiByNames(r,"CB","CG","CD","OE1"); break;
            case 'Q':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","CD");
                c3 = chiByNames(r,"CB","CG","CD","OE1"); break;
            case 'K':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","CD");
                c3 = chiByNames(r,"CB","CG","CD","CE");
                c4 = chiByNames(r,"CG","CD","CE","NZ"); break;
            case 'R':
                c1 = chiByNames(r,"N","CA","CB","CG");
                c2 = chiByNames(r,"CA","CB","CG","CD");
                c3 = chiByNames(r,"CB","CG","CD","NE");
                c4 = chiByNames(r,"CG","CD","NE","CZ"); break;
            default: break;             // G, A, X: no chi
        }
    }
};
