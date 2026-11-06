// NativeQCInputs.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Per-residue inputs for pdb_quality_check and structure_survey: the
// backbone dihedrals and the FoldX-style occupancy burial, computed from the
// coordinates.
// Built only from existing toolkit code:
//
//   phi / psi / omega
//     Residue::phi/psi/omega as set by Peptide::SetDihedrals(), which
//     PDB_Reader already runs on every chain at load. The -999 "unset"
//     sentinel becomes NaN. SetDihedrals detects only breaks that leave a gap
//     in the residue numbering (GAP residues have no atoms), so an angle that
//     spans a physical break (C(i-1)-N(i) > maxPeptideBond) is also set to NaN.
//
//   burial
//     FoldX's occupancy envelope: ComputeOccupancy() from FxAtoms.h, the
//     Gaussian-weighted volume of neighbouring heavy atoms (vol*exp(-d^2/24.5)
//     within faccCutoff, the residue's own atoms included), with FoldX per-atom
//     volumes, Occ and Occmax from FxLookup() (FoldX parameters_solvation
//     table). Every heavy atom of every chain is environment, so interface
//     residues count as buried. Residue burial is FoldX's ratio of sums over
//     the residue's atoms, clamped to [0,1]:
//         burial = sum(occupancy - Occ) / sum(Occmax - Occ)
//     scBurial over side-chain atoms, mcBurial over N, CA, C, O (OXT in
//     neither). A residue without side-chain atoms (Gly) gets scBurial = 0, as
//     FoldX writes it and as the reference DB bins expect (bin 0 holds all
//     glycines; a NaN would fall into FaccBin's last, fully buried, bin).
//     Each residue is scored against the atoms within (its radius + cutoff) of
//     its centroid, a strict superset of every atom's cutoff sphere, so the
//     result is identical to scoring against the whole structure.
//
//   Validated against FoldX 5 SequenceDetail on 12 structures (23,805
//   residues): dihedrals identical to 0.001 deg; occupancy rule identical,
//   values ~0.5% higher (the table is the FoldX 2014 parameter set).
//
// Nothing is dropped silently: NativeQCSummary counts skipped GAP residues,
// chain breaks, angles left undefined, and atoms that fell back to FxLookup's
// generic parameters (OXT counted separately: the FoldX table has no OXT).
#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "PDB_File.h"
#include "FxAtoms.h"       // FxAtom, FxAtomParam, FxLookup, ComputeOccupancy, Dist2

struct NativeResidueInputs
{
	int         pep = -1;          // index into pdb.peptides
	int         res = -1;          // index into pdb.peptides[pep].reslist
	char        chain = ' ';
	int         resnum = 0;        // Residue::number_in_pdb
	std::string resName;           // 3-letter, as read

	double phi = NAN, psi = NAN, omega = NAN;   // degrees; NaN = undefined

	double scBurial = NAN;         // side-chain burial in [0,1]; 0 when there are no side-chain atoms (Gly)
	double mcBurial = NAN;         // backbone (N, CA, C, O) burial in [0,1]; NaN only if none present
	int    nSc = 0, nMc = 0;       // atoms summed into scBurial / mcBurial
	int    nParamFallback = 0;     // this residue's atoms (OXT excluded) without a FoldX parameter

	std::vector<std::string> atomName;   // every heavy atom of the residue, in file order
	std::vector<double>      atomOcc;    // per-atom raw occupancy, parallel to atomName
};

struct NativeQCSummary
{
	int nResidues      = 0;   // residues returned
	int nGapResidues   = 0;   // GAP placeholders / atom-less residues skipped
	int nChainBreaks   = 0;   // consecutive residues with C(i-1)-N(i) > maxPeptideBond
	int nUnsetPhi      = 0;   // NaN after all rules (chain starts, gaps, breaks, missing atoms)
	int nUnsetPsi      = 0;
	int nUnsetOmega    = 0;
	int nNoSideChain   = 0;   // residues given scBurial = 0 because they have no side-chain atoms
	int nParamFallback = 0;   // atoms (OXT excluded) that got FxLookup's generic parameters
	int nOxtFallback   = 0;   // OXT atoms (always generic parameters)
};

inline bool NativeIsBackboneAtom(const std::string& n)
{
	return n == "N" || n == "CA" || n == "C" || n == "O" || n == "OXT";
}

// FoldX-repaired files name histidine by tautomer (H1S / H2S); the FoldX
// parameter table is keyed on HIS only.
inline const std::string& NativeParamResName(const std::string& r)
{
	static const std::string his = "HIS";
	return (r == "H1S" || r == "H2S") ? his : r;
}

inline std::vector<NativeResidueInputs>
ComputeNativeQCInputs(const PDB_File& pdb,
                      NativeQCSummary* summary = nullptr,
                      double faccCutoff = 6.0,       // FoldX occupancy cutoff (validated)
                      double maxPeptideBond = 2.0)   // A; ideal C-N is 1.33
{
	using ZipperScore::FxAtom;
	using ZipperScore::FxAtomParam;
	using ZipperScore::FxLookup;
	using ZipperScore::ComputeOccupancy;
	using ZipperScore::Dist2;

	NativeQCSummary s;
	// FxLookup returns one static object for every key it does not know.
	const FxAtomParam* fallback = &FxLookup("__none__", "__none__");

	// ── 1. Environment: every heavy atom of every chain, FoldX parameters ──
	std::vector<FxAtom> env;
	for (const auto& pep : pdb.peptides)
		for (const auto& r : pep.reslist)
		{
			if (r.residue_name == "GAP") continue;
			for (const auto& a : r.atoms)
			{
				FxAtom fa;
				fa.name  = a.atom_name;
				fa.pos   = Vector3D(a.x, a.y, a.z);
				fa.param = FxLookup(NativeParamResName(r.residue_name), a.atom_name);
				env.push_back(fa);
			}
		}

	// C(i-1)-N(i) peptide bond longer than maxPeptideBond = physical break.
	auto brokenBefore = [&](const std::vector<Residue>& rl, int i) -> bool
	{
		if (i <= 0) return false;
		const Residue& a = rl[i - 1];
		const Residue& b = rl[i];
		if (a.residue_name == "GAP" || b.residue_name == "GAP") return false;
		if (a.C < 0 || b.N < 0) return false;
		const Atom& c = a.atoms[a.C];
		const Atom& n = b.atoms[b.N];
		double dx = c.x - n.x, dy = c.y - n.y, dz = c.z - n.z;
		return dx * dx + dy * dy + dz * dz > maxPeptideBond * maxPeptideBond;
	};
	auto unset = [](double v) { return v <= -998.0 ? NAN : v; };   // SetDihedrals sentinel -999
	auto ratio = [](double occ, double lo, double hi) -> double     // FoldX ratio of sums, clamped
	{
		const double d = hi - lo;
		if (d <= 0.0) return NAN;
		return std::max(0.0, std::min(1.0, (occ - lo) / d));
	};

	std::vector<NativeResidueInputs> out;
	std::vector<FxAtom> local;

	for (int p = 0; p < (int)pdb.peptides.size(); ++p)
	{
		const std::vector<Residue>& rl = pdb.peptides[p].reslist;
		for (int i = 0; i < (int)rl.size(); ++i)
		{
			const Residue& r = rl[i];
			if (r.residue_name == "GAP" || r.atoms.empty()) { s.nGapResidues++; continue; }

			NativeResidueInputs q;
			q.pep = p; q.res = i;
			q.chain = r.chain; q.resnum = r.number_in_pdb; q.resName = r.residue_name;

			// ── 2. Dihedrals ──
			q.phi = unset(r.phi); q.psi = unset(r.psi); q.omega = unset(r.omega);
			if (brokenBefore(rl, i))      { q.phi = NAN; q.omega = NAN; s.nChainBreaks++; }
			if (i + 1 < (int)rl.size() && brokenBefore(rl, i + 1)) q.psi = NAN;
			if (std::isnan(q.phi))   s.nUnsetPhi++;
			if (std::isnan(q.psi))   s.nUnsetPsi++;
			if (std::isnan(q.omega)) s.nUnsetOmega++;

			// ── 3. Burial: local environment, occupancy per atom, ratio of sums ──
			double cx = 0, cy = 0, cz = 0;
			for (const auto& a : r.atoms) { cx += a.x; cy += a.y; cz += a.z; }
			cx /= r.atoms.size(); cy /= r.atoms.size(); cz /= r.atoms.size();
			const Vector3D centre(cx, cy, cz);
			double rmax2 = 0;
			for (const auto& a : r.atoms)
				rmax2 = std::max(rmax2, Dist2(centre, Vector3D(a.x, a.y, a.z)));
			const double R = std::sqrt(rmax2) + faccCutoff;
			local.clear();
			for (const auto& e : env)
				if (Dist2(centre, e.pos) <= R * R) local.push_back(e);

			double scOcc = 0, scLo = 0, scHi = 0, mcOcc = 0, mcLo = 0, mcHi = 0;
			for (const auto& a : r.atoms)
			{
				const FxAtomParam& prm = FxLookup(NativeParamResName(r.residue_name), a.atom_name);
				if (&prm == fallback)
				{
					if (a.atom_name == "OXT") s.nOxtFallback++;
					else { q.nParamFallback++; s.nParamFallback++; }
				}
				const double occ = ComputeOccupancy(Vector3D(a.x, a.y, a.z), local, faccCutoff);
				q.atomName.push_back(a.atom_name);
				q.atomOcc.push_back(occ);
				if (a.atom_name == "OXT") continue;   // terminal oxygen: in neither sum
				if (NativeIsBackboneAtom(a.atom_name)) { mcOcc += occ; mcLo += prm.Occ; mcHi += prm.Occmax; q.nMc++; }
				else                                    { scOcc += occ; scLo += prm.Occ; scHi += prm.Occmax; q.nSc++; }
			}
			if (q.nSc > 0) q.scBurial = ratio(scOcc, scLo, scHi);
			else         { q.scBurial = 0.0; s.nNoSideChain++; }
			if (q.nMc > 0) q.mcBurial = ratio(mcOcc, mcLo, mcHi);

			out.push_back(q);
			s.nResidues++;
		}
	}

	if (summary) *summary = s;
	return out;
}
