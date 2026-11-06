// Clashscore.cpp
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// See Clashscore.h for the description.
//
// Connectivity reuse: the covalent bond graph (and the 1-2/1-3/1-4 exclusion)
// uses ZipperScore::BuildResidueAdjacency, extended across the C(i)-N(i+1)
// peptide bond so
// junction pairs are excluded by true bond-path distance, not name heuristics.
// Hydrogens are added with ZipperScore::AddVirtualHydrogens.  Rotatable polar
// hydrogens (sidechain -OH / -SH / -NH3 / -NH2 / His), whose orientation is not
// determined without Reduce-style optimisation, are EXCLUDED from the clash
// count; backbone amide H and all C-H are kept.

#include "Clashscore.h"
#include "PDB_File.h"
#include "Peptide.h"
#include "Residue.h"
#include "Atom.h"
#include "vector3D.h"
#include "FxAtoms.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_set>

using namespace std;

namespace {

// Per-residue staging: heavy + virtual H atoms, their bond adjacency, and the
// mapping from local atom index → global clash-list index (-1 if the atom is
// excluded from the count, e.g. a rotatable polar H).
struct ResAtoms {
    vector<ZipperScore::FxAtom> fx;
    vector<vector<int>>         adj;
    vector<int>                 gidx;
    int  chainId;
    int  idxInChain;
    char chain;
    int  resnum;
    char aa1;
    string res3;
};

// One atom in the global clash list.
struct CAtom {
    float  x, y, z;
    float  r;        // clash VdW radius (min_radius)
    bool   polar;    // polar atom (for polar-polar tolerance)
    bool   polarH;   // backbone amide H (H-bond donor proton we keep)
    bool   acceptor; // heavy O / S, His ring N
    int    rgi;      // residue global index
    int    parent;   // for an H: global index of its heavy atom (-1 otherwise)
};

inline string trimName(const string& s) {
    size_t b = s.find_first_not_of(' ');
    if (b == string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

inline bool isHydrogenName(const string& n) {
    if (n.empty()) return false;
    if (n[0] == 'H') return true;
    if (n.size() >= 2 && isdigit((unsigned char)n[0]) && n[1] == 'H') return true;
    return false;
}

// probe radii: ClashRadii::ProbeRadius / isAromaticRingAtom (Clashscore.h), shared with the
// density-guided legality check so both use the same radii.
using ClashRadii::isAromaticRingAtom;
using ClashRadii::ProbeRadius;

// All-pairs bond-path distance within an adjacency list, capped at 4
// (we only care whether a pair is < 4 bonds apart).
inline void bfsDistLE4(const vector<vector<int>>& adj, vector<vector<int>>& dist) {
    int n = (int)adj.size();
    dist.assign(n, vector<int>(n, 99));
    for (int s = 0; s < n; ++s) {
        dist[s][s] = 0;
        vector<int> q{s}; size_t h = 0;
        while (h < q.size()) {
            int c = q[h++];
            if (dist[s][c] >= 4) continue;          // no need past 1-4
            for (int nb : adj[c]) if (dist[s][nb] == 99) {
                dist[s][nb] = dist[s][c] + 1;
                q.push_back(nb);
            }
        }
    }
}

// ── Movable-hydrogen (Reduce-style) optimisation ─────────────────────────────
// Sidechain amide / His-ring N-H are defined once the amide/ring flip state is
// chosen (they are NOT freely rotatable), so they are counted.
// Arg guanidinium NH: fixed by the planar heavy atoms, so counted (reduce builds them, probe
// counts them; they were dropped before 2026-10-07).
inline bool isPlanarArgNH(char aa, const string& nm) {
    return aa == 'R' && (nm == "HE" || nm == "HH11" || nm == "HH12" || nm == "HH21" || nm == "HH22");
}
inline bool isFlipDeterminedNH(char aa, const string& nm) {
    if (aa == 'N') return nm == "HD21" || nm == "HD22";
    if (aa == 'Q') return nm == "HE21" || nm == "HE22";
    if (aa == 'H') return nm == "HD1"  || nm == "HE2";
    return false;
}

// 180° rotation of p about the axis through a and b (the amide/ring flip).
inline Vector3D Flip180(const Vector3D& p, const Vector3D& a, const Vector3D& b) {
    double nx = b.x - a.x, ny = b.y - a.y, nz = b.z - a.z;
    double L = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (L < 1e-9) return p;
    nx /= L; ny /= L; nz /= L;
    double vx = p.x - a.x, vy = p.y - a.y, vz = p.z - a.z;
    double d = vx*nx + vy*ny + vz*nz;
    return Vector3D(2.0*(a.x + d*nx) - p.x,
                    2.0*(a.y + d*ny) - p.y,
                    2.0*(a.z + d*nz) - p.z);
}

// Rodrigues rotation of p by angle (cos=c, sin=s) about the unit axis (nx,ny,nz)
// through the point o.
inline Vector3D RotateAxis(const Vector3D& p, const Vector3D& o,
                           double nx, double ny, double nz, double c, double s) {
    double vx = p.x - o.x, vy = p.y - o.y, vz = p.z - o.z;
    double dot = vx*nx + vy*ny + vz*nz;
    double cx = ny*vz - nz*vy, cy = nz*vx - nx*vz, cz = nx*vy - ny*vx; // n × v
    return Vector3D(o.x + vx*c + cx*s + nx*dot*(1.0 - c),
                    o.y + vy*c + cy*s + ny*dot*(1.0 - c),
                    o.z + vz*c + cz*s + nz*dot*(1.0 - c));
}

// Spin rotor: H on `center`, spun about the center→heavy bond.  3-fold groups
// (methyl, -NH3) search 0..120°; single H (-OH, -SH) search 0..360°.
struct SpinRotor { const char* center; const char* heavy; std::vector<std::string> h; int searchDeg; };
inline std::vector<SpinRotor> spinRotors(char aa) {
    // Carbon methyls (Ala, Val, Thr CG2, Ile, Leu) are NOT spun: reduce -build leaves them
    // staggered (checked on its output: every Ala/Val methyl H at exactly 180/-60/60), so
    // spinning them hid clashes probe counts.  Met CE kept (not checked: no Met tested).
    switch (aa) {
        case 'T': return {{"OG1","CB",{"HG1"},360}};                 // hydroxyl
        case 'M': return {{"CE","SD",{"HE1","HE2","HE3"},120}};
        case 'S': return {{"OG","CB",{"HG"},360}};                   // hydroxyl
        case 'Y': return {{"OH","CZ",{"HH"},360}};                   // hydroxyl
        case 'C': return {{"SG","CB",{"HG"},360}};                   // thiol
        case 'K': return {{"NZ","CE",{"HZ1","HZ2","HZ3"},120}};      // ammonium
        default:  return {};
    }
}

// Rotatable polar sidechain H that this optimiser places (so they are counted,
// not dropped).  Arg guanidinium H are not yet handled and stay dropped.
inline bool isOptimizedPolarH(char aa, const std::string& nm) {
    if (aa == 'S') return nm == "HG";
    if (aa == 'T') return nm == "HG1";
    if (aa == 'Y') return nm == "HH";
    if (aa == 'C') return nm == "HG";
    if (aa == 'K') return nm == "HZ1" || nm == "HZ2" || nm == "HZ3";
    return false;
}

// Set each movable-H rotor to its lowest-clash orientation, as Reduce does
// before MolProbity counts.  Currently: Asn/Gln amide flip and His ring flip
// (methyls and -OH/-SH/-NH3 rotors land in follow-up commits).  Judged by the
// inter-residue clash-overlap the rotor atoms incur; iterated so flips on
// stacked copies converge.  Operates in place on R.fx positions.
inline void OptimizeMovableH(std::vector<ResAtoms>& R, double overlapCut) {
    struct FA { float x, y, z, r; int res; };
    const int ITER = 3;
    // Flat-array base offset of each residue (atom counts are constant across
    // iterations; only positions change), so a residue-local index maps to its
    // flat index as resBase[r] + localIndex.
    vector<int> resBase(R.size() + 1, 0);
    for (int r = 0; r < (int)R.size(); ++r)
        resBase[r+1] = resBase[r] + (int)R[r].fx.size();

    for (int it = 0; it < ITER; ++it) {
        vector<FA> fa; fa.reserve(resBase.back());
        for (int r = 0; r < (int)R.size(); ++r)
            for (const auto& f : R[r].fx)
                fa.push_back(FA{(float)f.pos.x,(float)f.pos.y,(float)f.pos.z,
                                (float)f.param.min_radius, r});
        if (fa.size() < 2) return;

        float xmin=fa[0].x,xmax=xmin,ymin=fa[0].y,ymax=ymin,zmin=fa[0].z,zmax=zmin;
        for (const auto& a : fa) { xmin=min(xmin,a.x);xmax=max(xmax,a.x);
            ymin=min(ymin,a.y);ymax=max(ymax,a.y);zmin=min(zmin,a.z);zmax=max(zmax,a.z); }
        const float cell = 4.0f;
        int nx=max(1,(int)((xmax-xmin)/cell)+1), ny=max(1,(int)((ymax-ymin)/cell)+1),
            nz=max(1,(int)((zmax-zmin)/cell)+1);
        auto cellOf=[&](float x,float y,float z,int&ix,int&iy,int&iz){
            ix=(int)((x-xmin)/cell); if(ix<0)ix=0; if(ix>=nx)ix=nx-1;
            iy=(int)((y-ymin)/cell); if(iy<0)iy=0; if(iy>=ny)iy=ny-1;
            iz=(int)((z-zmin)/cell); if(iz<0)iz=0; if(iz>=nz)iz=nz-1; };
        vector<vector<int>> bins((size_t)nx*ny*nz);
        for (int i=0;i<(int)fa.size();++i){ int ix,iy,iz;
            cellOf(fa[i].x,fa[i].y,fa[i].z,ix,iy,iz);
            bins[(size_t)ix+(size_t)nx*(iy+(size_t)ny*iz)].push_back(i); }

        // Clash-overlap of a candidate atom (pos,radius) vs all atoms except
        // those in `skip` (flat indices bonded within 4 of the rotor, incl. the
        // rotor's own atoms).  Same-residue siblings (e.g. the other methyl of
        // a Val) are NOT skipped, so the optimiser cannot spin a rotor into its
        // own residue.
        auto atomClash=[&](double px,double py,double pz,double rad,
                           const std::unordered_set<int>& skip)->double{
            int ix,iy,iz; cellOf((float)px,(float)py,(float)pz,ix,iy,iz); double s=0;
            for(int dz=-1;dz<=1;++dz)for(int dy=-1;dy<=1;++dy)for(int dx=-1;dx<=1;++dx){
                int jx=ix+dx,jy=iy+dy,jz=iz+dz;
                if(jx<0||jy<0||jz<0||jx>=nx||jy>=ny||jz>=nz)continue;
                for(int j:bins[(size_t)jx+(size_t)nx*(jy+(size_t)ny*jz)]){
                    if(skip.count(j))continue;
                    double ex=px-fa[j].x,ey=py-fa[j].y,ez=pz-fa[j].z;
                    double d2=ex*ex+ey*ey+ez*ez, rsum=rad+fa[j].r;
                    if(d2>=rsum*rsum)continue;
                    double ov=rsum-std::sqrt(d2); if(ov>overlapCut)s+=ov;
                }
            }
            return s;
        };
        // Flat indices within residue r bonded within 4 of any moving atom.
        auto buildSkip=[&](int r,const vector<vector<int>>& dist,
                           const vector<int>& moving)->std::unordered_set<int>{
            std::unordered_set<int> skip; int m=(int)R[r].fx.size();
            for(int b=0;b<m;++b)
                for(int a:moving) if(dist[a][b]<4){ skip.insert(resBase[r]+b); break; }
            return skip;
        };

        bool anyChange=false; int nFlip=0; double sumCur=0,sumFlp=0;
        for (int r=0;r<(int)R.size();++r){
            char aa=R[r].aa1;
            if(aa!='N'&&aa!='Q'&&aa!='H')continue;
            ResAtoms& ra=R[r];
            const char* axA; const char* axB; vector<string> grp;
            if(aa=='N'){axA="CB";axB="CG";grp={"OD1","ND2","HD21","HD22"};}
            else if(aa=='Q'){axA="CG";axB="CD";grp={"OE1","NE2","HE21","HE22"};}
            else {axA="CB";axB="CG";grp={"ND1","CD2","CE1","NE2","HD1","HD2","HE1","HE2"};}
            auto findPos=[&](const char* nm,Vector3D& o)->bool{
                for(const auto&f:ra.fx)if(f.name==nm){o=f.pos;return true;} return false; };
            Vector3D A,B; if(!findPos(axA,A)||!findPos(axB,B))continue;
            vector<int> gi;
            for(int li=0;li<(int)ra.fx.size();++li)
                for(const auto&gn:grp)if(ra.fx[li].name==gn){gi.push_back(li);break;}
            if(gi.empty())continue;
            vector<vector<int>> dist; bfsDistLE4(ra.adj,dist);
            auto skip=buildSkip(r,dist,gi);
            double cur=0; for(int li:gi)
                cur+=atomClash(ra.fx[li].pos.x,ra.fx[li].pos.y,ra.fx[li].pos.z,
                               ra.fx[li].param.min_radius,skip);
            double flp=0; vector<Vector3D> np(gi.size());
            for(size_t k=0;k<gi.size();++k){
                np[k]=Flip180(ra.fx[gi[k]].pos,A,B);
                flp+=atomClash(np[k].x,np[k].y,np[k].z,ra.fx[gi[k]].param.min_radius,skip);
            }
            if(flp<cur-1e-3){ for(size_t k=0;k<gi.size();++k)ra.fx[gi[k]].pos=np[k];
                              anyChange=true; ++nFlip; sumCur+=cur; sumFlp+=flp; }
        }

        // ── spin rotors: methyl / -NH3 / -OH / -SH, about the center–heavy bond ──
        int nSpin=0;
        for (int r=0;r<(int)R.size();++r){
            ResAtoms& ra=R[r];
            auto rotors=spinRotors(ra.aa1);
            if(rotors.empty())continue;
            vector<vector<int>> dist; bfsDistLE4(ra.adj,dist);
            auto findPos=[&](const char* nm,Vector3D& o)->bool{
                for(const auto&f:ra.fx)if(f.name==nm){o=f.pos;return true;} return false; };
            for(const auto& m : rotors){
                Vector3D C,Hv; if(!findPos(m.center,C)||!findPos(m.heavy,Hv))continue;
                double nx=C.x-Hv.x,ny=C.y-Hv.y,nz=C.z-Hv.z;
                double L=std::sqrt(nx*nx+ny*ny+nz*nz); if(L<1e-9)continue;
                nx/=L; ny/=L; nz/=L;
                vector<int> hi;
                for(const auto& hn : m.h)
                    for(int li=0;li<(int)ra.fx.size();++li)
                        if(ra.fx[li].name==hn){ hi.push_back(li); break; }
                if(hi.size()!=m.h.size())continue;
                auto skip=buildSkip(r,dist,hi);
                double bestDelta=0, bestClash=1e18;
                for(int deg=0; deg<m.searchDeg; deg+=10){
                    double th=deg*M_PI/180.0, c=std::cos(th), s=std::sin(th), sum=0;
                    for(int li:hi){
                        Vector3D q=RotateAxis(ra.fx[li].pos,C,nx,ny,nz,c,s);
                        sum+=atomClash(q.x,q.y,q.z,ra.fx[li].param.min_radius,skip);
                    }
                    if(sum<bestClash-1e-6){bestClash=sum; bestDelta=th;}
                }
                if(bestDelta!=0){
                    double c=std::cos(bestDelta), s=std::sin(bestDelta);
                    for(int li:hi)
                        ra.fx[li].pos=RotateAxis(ra.fx[li].pos,C,nx,ny,nz,c,s);
                    ++nSpin;
                }
            }
        }

        if(getenv("CLASH_DBG"))
            fprintf(stderr,"[movableH] iter %d: %d flips (overlap %.1f->%.1f), %d methyl spins\n",
                    it,nFlip,sumCur,sumFlp,nSpin);
        if(!anyChange && nSpin==0)break;
    }
}

} // namespace

std::vector<ClashResult> Clashscore::Compute(const PDB_File& pdb,
                                             double overlapCut,
                                             double hbondAllow,
                                             bool   includeHydrogens)
{
    vector<ClashResult> out;
    vector<ResAtoms>    R;
    vector<CAtom>       atoms;
    vector<string>      lbl;          // parallel to atoms, for CLASH_PAIRS dump
    const bool dumpPairs = getenv("CLASH_PAIRS") != nullptr;   // dev diagnostic
    atoms.reserve(8192);

    int rgi = 0;
    for (size_t pi = 0; pi < pdb.peptides.size(); ++pi)
    {
        const Peptide& pep = pdb.peptides[pi];
        int idxInChain = 0;
        for (const Residue& res : pep.reslist)
        {
            if (res.residue_name == "HOH" || res.residue_name == "GAP") continue;
            char aa1 = res.OneLetter();

            ResAtoms ra;
            ra.chainId = (int)pi; ra.idxInChain = idxInChain;
            ra.chain = res.chain; ra.resnum = res.number_in_pdb; ra.aa1 = aa1;
            ra.res3 = res.residue_name;

            // Heavy atoms (strip any deposited H for consistency).
            for (const Atom& a : res.atoms) {
                string nm = trimName(a.atom_name);
                if (nm.empty() || isHydrogenName(nm)) continue;
                if (nm.compare(0, 2, "X_") == 0) continue;   // pseudo-atom (aromatic centroid), not an atom
                ZipperScore::FxAtom f;
                f.name = nm;
                f.pos  = Vector3D(a.x, a.y, a.z);
                f.param = ZipperScore::FxAtomParam{};
                string elem = ZipperScore::ElementFromName(nm);
                f.param.min_radius = ZipperScore::VdwRadius(elem);
                f.param.radius     = f.param.min_radius;
                f.param.polar      = (elem == "N" || elem == "O" || elem == "S") ? 1 : 0;
                f.resSeqIdx = rgi;
                ra.fx.push_back(f);
            }
            if (ra.fx.empty()) continue;

            if (includeHydrogens && aa1 != 'X')
                ZipperScore::AddVirtualHydrogens(ra.fx, aa1, res.residue_name);

            // Override the force-field VdW radii with MolProbity/probe radii,
            // locally, for the clash metric only (heavy atoms re-set from line
            // 128; H radii set here since AddVirtualHydrogens uses 0.99).
            for (auto& f : ra.fx) {
                bool isH    = isHydrogenName(f.name);
                string elem = isH ? "H" : ZipperScore::ElementFromName(f.name);
                bool polarH = isH && (f.param.polar == 1);
                f.param.min_radius = ProbeRadius(aa1, f.name, elem, polarH);
                f.param.radius     = f.param.min_radius;
            }

            // Bond adjacency (shared with TorsionalClash).
            if (!ZipperScore::BuildResidueAdjacency(ra.fx, aa1, ra.adj))
                ra.adj.assign(ra.fx.size(), {});   // unknown AA → no intra exclusions

            // Correct H→parent bonds geometrically (each H bonds to its nearest
            // heavy atom).  BuildResidueAdjacency's name heuristic mis-bonds the
            // sidechain amide / His-ring NH (e.g. Asn HD21 → backbone C instead
            // of ND2, since Asn has no CD2), which would leave the N-H bond
            // itself counted as a 1.0 Å clash.  Local to this metric; the shared
            // adjacency the force field uses is untouched.
            if (includeHydrogens) {
                for (int li = 0; li < (int)ra.fx.size(); ++li) {
                    if (!isHydrogenName(ra.fx[li].name)) continue;
                    int par = -1; double best = 1e18;
                    for (int hj = 0; hj < (int)ra.fx.size(); ++hj) {
                        if (hj == li || isHydrogenName(ra.fx[hj].name)) continue;
                        double dx = ra.fx[li].pos.x - ra.fx[hj].pos.x;
                        double dy = ra.fx[li].pos.y - ra.fx[hj].pos.y;
                        double dz = ra.fx[li].pos.z - ra.fx[hj].pos.z;
                        double d2 = dx*dx + dy*dy + dz*dz;
                        if (d2 < best) { best = d2; par = hj; }
                    }
                    if (par < 0) continue;
                    // reduce -build bond lengths (measured on its output, 2026-10-07):
                    // C-H 0.97, aromatic C-H 0.93, N-H 0.86, NH3+ 0.89, O-H 0.84
                    {   const string& pn = ra.fx[par].name; string pe = ZipperScore::ElementFromName(pn);
                        double L = 0;
                        if (pe == "C") L = isAromaticRingAtom(aa1, pn) ? 0.93 : 0.97;
                        else if (pe == "N") L = (aa1 == 'K' && pn == "NZ") ? 0.89 : 0.86;
                        else if (pe == "O") L = 0.84;
                        if (L > 0) { Vector3D& hp = ra.fx[li].pos; const Vector3D& pp = ra.fx[par].pos;
                            Vector3D d(hp.x - pp.x, hp.y - pp.y, hp.z - pp.z);
                            double n = std::sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
                            if (n > 1e-6) hp = Vector3D(pp.x + d.x * L / n, pp.y + d.y * L / n, pp.z + d.z * L / n); } }
                    for (auto& lst : ra.adj)
                        lst.erase(std::remove(lst.begin(), lst.end(), li), lst.end());
                    ra.adj[li].clear();
                    ra.adj[li].push_back(par);
                    ra.adj[par].push_back(li);
                }
            }

            // CAtom construction is deferred to a post-loop pass, so the
            // movable-H optimiser (OptimizeMovableH) can first set each rotor
            // to its lowest-clash orientation.  R index == rgi (both advance in
            // lockstep over accepted residues).
            out.push_back(ClashResult{res.chain, res.number_in_pdb, aa1, 0.0});
            R.push_back(std::move(ra));
            ++rgi; ++idxInChain;
        }
    }

    // ── reduce-style movable-hydrogen optimisation (with-H path only) ─────────
    // Spin methyls / -OH / -SH / -NH3 and flip Asn/Gln/His to the orientation
    // with the fewest clashes BEFORE counting, exactly as MolProbity's Reduce
    // does.  Counting fixed-idealised rotatable H otherwise inflates the score
    // with clashes a simple rotation removes (and the optimiser then chases
    // those ghosts).  Local to the metric; the force field is untouched.
    // ── hydrogens as reduce -build places them (metric-local; force field untouched) ──
    // (1) backbone amide H on the bisector of C(i-1)-N and CA-N at 0.86 A.  AddVirtualHydrogens
    //     only sees one residue and builds it from the residue's OWN C, ~1.0 A off reduce's H.
    // (2) carbon methyls staggered as reduce leaves them: H1 anti (180) to the reference atom,
    //     H2 at -60, H3 at +60 (checked on reduce output: Ala N-CA-CB-HB1/2/3 = 180/-60/60,
    //     Val CA-CB-CG1-HG11 = 180), C-H 0.97, angle 109.5.
    if (includeHydrogens) {
        auto idx = [](const ResAtoms& x, const char* nm) {
            for (int k = 0; k < (int)x.fx.size(); ++k) if (x.fx[k].name == nm) return k; return -1; };
        for (size_t r = 0; r < R.size(); ++r) {
            ResAtoms& ra = R[r];
            int iH = idx(ra, "H"), iN = idx(ra, "N"), iCA = idx(ra, "CA");
            if (iH >= 0 && iN >= 0 && iCA >= 0 && r > 0 && R[r-1].chainId == ra.chainId &&
                R[r-1].idxInChain + 1 == ra.idxInChain) {
                int iC = idx(R[r-1], "C");
                if (iC >= 0) {
                    const Vector3D& N = ra.fx[iN].pos; const Vector3D& CA = ra.fx[iCA].pos; const Vector3D& Cp = R[r-1].fx[iC].pos;
                    Vector3D u1 = N - Cp, u2 = N - CA; u1.Normalize(); u2.Normalize();
                    Vector3D h = u1 + u2;
                    if (std::sqrt(h * h) > 1e-6 && std::sqrt((N - Cp) * (N - Cp)) < 2.0) {
                        h.Normalize(); ra.fx[iH].pos = N + 0.86 * h; }
                }
            }
            struct Me { const char* c; const char* x; const char* ref; const char* h[3]; };
            static const std::vector<std::pair<char, std::vector<Me>>> methyls = {
                {'A', {{"CB","CA","N",{"HB1","HB2","HB3"}}}},
                {'V', {{"CG1","CB","CA",{"HG11","HG12","HG13"}}, {"CG2","CB","CA",{"HG21","HG22","HG23"}}}},
                {'T', {{"CG2","CB","CA",{"HG21","HG22","HG23"}}}},
                {'I', {{"CG2","CB","CA",{"HG21","HG22","HG23"}}, {"CD1","CG1","CB",{"HD11","HD12","HD13"}}}},
                {'L', {{"CD1","CG","CB",{"HD11","HD12","HD13"}}, {"CD2","CG","CB",{"HD21","HD22","HD23"}}}}};
            for (const auto& m : methyls) {
                if (m.first != ra.aa1) continue;
                for (const auto& g : m.second) {
                    int ic = idx(ra, g.c), ix = idx(ra, g.x), ir = idx(ra, g.ref);
                    if (ic < 0 || ix < 0 || ir < 0) continue;
                    static const double iupac[3] = {180.0, -60.0, 60.0};
                    for (int k = 0; k < 3; ++k) {
                        int ih = idx(ra, g.h[k]); if (ih < 0) continue;
                        ra.fx[ih].pos = ZipperScore::PlaceAtom(ra.fx[ir].pos, ra.fx[ix].pos, ra.fx[ic].pos,
                                                               0.97, 109.5, 180.0 - iupac[k]);  // PlaceAtom: 180 - IUPAC
                    }
                }
            }
        }
    }

    if (includeHydrogens)
        OptimizeMovableH(R, overlapCut);

    // Dev diagnostic: dump the exact atom set we count (heavy + kept H, after
    // movable-H optimisation) so probe can be run on our hydrogenated model.
    if (includeHydrogens && getenv("CLASH_DUMP")) {
        FILE* fp = fopen(getenv("CLASH_DUMP"), "w");
        if (fp) {
            int serial = 1;
            for (const ResAtoms& ra : R) {
                for (const auto& f : ra.fx) {
                    bool isH    = isHydrogenName(f.name);
                    bool polar  = (f.param.polar == 1);
                    bool amideH = (f.name == "H" || f.name == "HN");
                    bool flipNH  = isFlipDeterminedNH(ra.aa1, f.name) || isPlanarArgNH(ra.aa1, f.name);
                    bool optPolH = isOptimizedPolarH(ra.aa1, f.name);
                    if (isH && polar && !amideH && !flipNH && !optPolH) continue;  // same keep-rule
                    string nm = f.name;
                    string nmField = (nm.size() >= 4) ? nm : " " + nm;
                    while (nmField.size() < 4) nmField += " ";
                    fprintf(fp, "ATOM  %5d %-4s %-3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                            serial++, nmField.c_str(), ra.res3.c_str(), ra.chain,
                            ra.resnum, f.pos.x, f.pos.y, f.pos.z);
                }
            }
            fclose(fp);
        }
    }

    // ── Build the global clash list from the (now optimised) atoms ────────────
    for (int rr = 0; rr < (int)R.size(); ++rr) {
        ResAtoms& ra = R[rr];
        ra.gidx.assign(ra.fx.size(), -1);
        for (int li = 0; li < (int)ra.fx.size(); ++li) {
            const auto& f = ra.fx[li];
            bool isH    = isHydrogenName(f.name);
            bool polar  = (f.param.polar == 1);
            bool amideH = (f.name == "H" || f.name == "HN");
            // Flip-determined polar H (sidechain amide NH2 / His ring NH) are
            // NOT freely rotatable: once the flip state is chosen they are
            // defined, so they ARE counted.  Truly rotatable polar H (-OH /
            // -SH / -NH3) remain skipped until their rotor optimiser lands.
            bool flipNH  = isFlipDeterminedNH(ra.aa1, f.name) || isPlanarArgNH(ra.aa1, f.name);
            bool optPolH = isOptimizedPolarH(ra.aa1, f.name);
            // Drop only polar H we have NOT placed definitively: Arg guanidinium
            // etc.  Backbone amide, flip-determined NH, and the optimised
            // rotatable -OH/-SH/-NH3 are all counted (matching probe).
            if (isH && polar && !amideH && !flipNH && !optPolH) continue;
            // (heavy-only mode never adds H, so isH is false there)
            CAtom c;
            c.x = (float)f.pos.x; c.y = (float)f.pos.y; c.z = (float)f.pos.z;
            c.r = (float)f.param.min_radius;
            c.polar = polar;
            c.polarH = isH && (amideH || flipNH || optPolH);  // H-bond donor protons
            c.acceptor = (!isH) &&
                [&]{ string e = ZipperScore::ElementFromName(f.name);
                     if (e == "O" || e == "S") return true;
                     return ra.aa1 == 'H' && (f.name == "ND1" || f.name == "NE2"); }();   // not amide/amine N
            c.rgi = rr;
            c.parent = -1;
            ra.gidx[li] = (int)atoms.size();
            atoms.push_back(c);
            if (dumpPairs)
                lbl.push_back(string(1,ra.chain)+std::to_string(ra.resnum)+ra.res3+"/"+f.name);
        }
    }

    // H -> heavy-atom parent (global indices), from the residue adjacency
    for (const auto& ra : R)
        for (int li = 0; li < (int)ra.fx.size(); ++li) {
            int gi = ra.gidx[li];
            if (gi < 0 || !isHydrogenName(ra.fx[li].name) || ra.adj[li].empty()) continue;
            atoms[gi].parent = ra.gidx[ra.adj[li][0]];
        }
    const int nAtoms = (int)atoms.size();
    if (nAtoms < 2) return out;

    // ── Build excluded-pair set ──────────────────────────────────────────────
    // probe (checked on its own output, 2026-10-07): pairs <= 3 bonds apart are never
    // reported; 4-bond pairs are reported between heavy atoms (e.g. O(i)-C(i+1),
    // C-CD1 of Ile) but never when a hydrogen is involved.
    unordered_set<long long> excl;
    const long long N = nAtoms;
    auto addExcl = [&](int gi, int gj) {
        if (gi < 0 || gj < 0 || gi == gj) return;
        if (gi > gj) std::swap(gi, gj);
        excl.insert((long long)gi * N + gj);
    };

    // intra-residue
    for (const auto& ra : R) {
        vector<vector<int>> dist;
        bfsDistLE4(ra.adj, dist);
        int m = (int)ra.fx.size();
        for (int a = 0; a < m; ++a)
            for (int b = a + 1; b < m; ++b)
                if (dist[a][b] <= 3 || (dist[a][b] == 4 &&
                    (isHydrogenName(ra.fx[a].name) || isHydrogenName(ra.fx[b].name))))
                    addExcl(ra.gidx[a], ra.gidx[b]);   // probe, checked on its own output
    }

    // peptide-junction (consecutive residues, same chain): combined BFS with the
    // C(i)-N(i+1) bond added, so 1-2/1-3/1-4 pairs across the junction are excluded.
    for (size_t r = 0; r + 1 < R.size(); ++r) {
        const ResAtoms& A = R[r];
        const ResAtoms& B = R[r + 1];
        if (A.chainId != B.chainId || B.idxInChain != A.idxInChain + 1) continue;
        int na = (int)A.fx.size(), nb = (int)B.fx.size();
        // combined adjacency
        vector<vector<int>> cadj(na + nb);
        for (int i = 0; i < na; ++i) for (int j : A.adj[i]) cadj[i].push_back(j);
        for (int i = 0; i < nb; ++i) for (int j : B.adj[i]) cadj[na + i].push_back(na + j);
        auto find = [](const ResAtoms& x, const char* nm) {
            for (int i = 0; i < (int)x.fx.size(); ++i) if (x.fx[i].name == nm) return i;
            return -1;
        };
        int ci = find(A, "C"), nj = find(B, "N");
        if (ci < 0 || nj < 0) continue;            // chain break → no peptide bond
        cadj[ci].push_back(na + nj);
        cadj[na + nj].push_back(ci);
        vector<vector<int>> dist;
        bfsDistLE4(cadj, dist);
        for (int a = 0; a < na; ++a)
            for (int b = 0; b < nb; ++b)
                if (dist[a][na + b] <= 3 || (dist[a][na + b] == 4 &&
                    (isHydrogenName(A.fx[a].name) || isHydrogenName(B.fx[b].name))))
                    addExcl(A.gidx[a], B.gidx[b]);     // probe, checked on its own output
    }

    // ── Uniform spatial grid ────────────────────────────────────────────────
    float xmin = atoms[0].x, xmax = xmin, ymin = atoms[0].y, ymax = ymin,
          zmin = atoms[0].z, zmax = zmin;
    for (const auto& a : atoms) {
        xmin = min(xmin, a.x); xmax = max(xmax, a.x);
        ymin = min(ymin, a.y); ymax = max(ymax, a.y);
        zmin = min(zmin, a.z); zmax = max(zmax, a.z);
    }
    const float cell = 4.0f;
    int nx = max(1, (int)((xmax - xmin) / cell) + 1);
    int ny = max(1, (int)((ymax - ymin) / cell) + 1);
    int nz = max(1, (int)((zmax - zmin) / cell) + 1);
    auto cellOf = [&](float x, float y, float z, int& ix, int& iy, int& iz) {
        ix = (int)((x - xmin) / cell); if (ix < 0) ix = 0; if (ix >= nx) ix = nx - 1;
        iy = (int)((y - ymin) / cell); if (iy < 0) iy = 0; if (iy >= ny) iy = ny - 1;
        iz = (int)((z - zmin) / cell); if (iz < 0) iz = 0; if (iz >= nz) iz = nz - 1;
    };
    vector<vector<int>> bins((size_t)nx * ny * nz);
    for (int i = 0; i < nAtoms; ++i) {
        int ix, iy, iz; cellOf(atoms[i].x, atoms[i].y, atoms[i].z, ix, iy, iz);
        bins[(size_t)ix + (size_t)nx * (iy + (size_t)ny * iz)].push_back(i);
    }

    // ── Count clashing contacts ─────────────────────────────────────────────
    struct Hit { int i, j; double d, rsum, ov; };
    vector<Hit> found;
    for (int i = 0; i < nAtoms; ++i) {
        const CAtom& ai = atoms[i];
        int ix, iy, iz; cellOf(ai.x, ai.y, ai.z, ix, iy, iz);
        for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
            int jx = ix + dx, jy = iy + dy, jz = iz + dz;
            if (jx < 0 || jy < 0 || jz < 0 || jx >= nx || jy >= ny || jz >= nz) continue;
            for (int j : bins[(size_t)jx + (size_t)nx * (jy + (size_t)ny * jz)]) {
                if (j <= i) continue;
                const CAtom& aj = atoms[j];
                if (ai.rgi == aj.rgi) {
                    if (excl.count((long long)i * N + j)) continue;
                } else if (excl.count((long long)i * N + j)) {
                    continue;
                }
                float ex = ai.x - aj.x, ey = ai.y - aj.y, ez = ai.z - aj.z;
                float d2 = ex*ex + ey*ey + ez*ez;
                float rsum = ai.r + aj.r;
                if (d2 >= rsum * rsum) continue;
                double d = std::sqrt((double)d2);
                double overlap = (double)rsum - d;
                // NOTE: no polar-polar tolerance here: MolProbity/probe has no
                // such term; it only allows donor-H···acceptor H-bond overlap.
                // donor H / acceptor: an H-bond up to an overlap of hbondAllow (0.6), a
                // clash beyond it (probe: 'hb' never below a gap of -0.58; -0.65 is a clash)
                if ((ai.polarH && aj.acceptor) || (aj.polarH && ai.acceptor)) {
                    if (overlap <= hbondAllow) continue;
                }
                // disulfide bond (Cys SG-SG): bonded, not a clash
                if (overlap >= overlapCut) {
                    // Skip a Cys SG–SG disulfide (a bond, not a clash).  The
                    // test must identify SULFUR (radius 1.80); the old `> 1.7`
                    // also matched carbon (1.70f > 1.7 is true in float), which
                    // silently discarded every C–C contact under 2.5 Å, i.e.
                    // the most severe clashes in a tightly packed core.
                    if (d < 2.5 && ai.r > 1.75 && aj.r > 1.75) continue; // S-S bond
                    found.push_back({i, j, d, rsum, overlap});
                }
            }
        }
    }

    // NOTE: probe reports some heavy-atom contacts through a covering H instead (dot surface);
    // no simple pairwise rule reproduces that (tried 2026-10-07: blanket and geometric versions
    // both lose real pairs), so those contacts are counted on both pairs here (~3% of pairs).
    for (const auto& h : found) {
        out[atoms[h.i].rgi].clashscore += 1.0;
        out[atoms[h.j].rgi].clashscore += 1.0;
        if (dumpPairs && includeHydrogens)
            fprintf(stderr, "PAIR %-18s %-18s d=%.2f rsum=%.2f ov=%.2f%s\n",
                    lbl[h.i].c_str(), lbl[h.j].c_str(), h.d, h.rsum, h.ov,
                    atoms[h.i].rgi == atoms[h.j].rgi ? " INTRA" : "");
    }
    return out;
}

double Clashscore::StructureClashscore(const std::vector<ClashResult>& res, int nAtoms)
{
    if (nAtoms <= 0) return 0.0;
    double total = 0.0;
    for (const auto& r : res) total += r.clashscore;
    total *= 0.5;
    return total / nAtoms * 1000.0;
}
