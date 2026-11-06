// structure_survey_command.cpp
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz

#include "structure_survey_command.h"
#include "NativeQCInputs.h"
#include "BondGeometry.h"
#include "Clashscore.h"
#include "SurfacePatcher.h"
#include "CavityDetector.h"
#include "StructureReferenceDB.h"
#include "PDB_File.h"
#include "PDB_Reader.h"
#include "LineReader.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <map>

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
// Per-residue input record. Filled from NativeQCInputs.h (dihedrals and
// FoldX-style side-chain burial computed from the coordinates, exactly as in
// pdb_quality_check, so reference and query residues are binned identically).
// The energy metrics stay NaN, so the reference cells hold no energy
// distributions.
// ─────────────────────────────────────────────────────────────────────────────

struct SDParsed
{
    string residueName;  // 3-letter
    char   chain;
    int    resNumber;
    char   ss;           // '-' (not assigned)
    double facc;         // side-chain burial in [0,1], HIGH = buried
    double values[SM_COUNT];  // NaN-initialised

    bool valid = false;
};

// Convert 3-letter residue name to 1-letter
static char threeToOne(const string& three)
{
    if (three == "ALA") return 'A'; if (three == "CYS") return 'C';
    if (three == "ASP") return 'D'; if (three == "GLU") return 'E';
    if (three == "PHE") return 'F'; if (three == "GLY") return 'G';
    // FoldX writes histidine as HIS / H1S / H2S depending on the tautomer it
    // assigns; all three are histidine.
    if (three == "HIS" || three == "H1S" || three == "H2S") return 'H';
    if (three == "ILE") return 'I';
    if (three == "LYS") return 'K'; if (three == "LEU") return 'L';
    if (three == "MET") return 'M'; if (three == "ASN") return 'N';
    if (three == "PRO") return 'P'; if (three == "GLN") return 'Q';
    if (three == "ARG") return 'R'; if (three == "SER") return 'S';
    if (three == "THR") return 'T'; if (three == "VAL") return 'V';
    if (three == "TRP") return 'W'; if (three == "TYR") return 'Y';
    return '?';
}

// ═════════════════════════════════════════════════════════════════════════════
// Execute
// ═════════════════════════════════════════════════════════════════════════════
void structure_survey_command::Execute()
{
    if (Option.inputfiles.empty())
    {
        cerr << "structure_survey: no input file specified.\n";
        return;
    }

    // ── Read PDB list ────────────────────────────────────────────────────
    vector<string> pdbPaths;
    {
        ifstream list(Option.inputfiles[0].c_str());
        if (!list) {
            cerr << "structure_survey: cannot open " << Option.inputfiles[0] << endl;
            return;
        }
        string line;
        while (getline(list, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            // Prepend pdb dir if needed
            if (!Option.surface_pdb_dir.empty() && line[0] != '/' && line[0] != '~') {
                string dir = Option.surface_pdb_dir;
                if (dir.back() != '/') dir += '/';
                line = dir + line;
            }
            pdbPaths.push_back(line);
        }
    }

    int nPdbs = (int)pdbPaths.size();
    cerr << "structure_survey: " << nPdbs << " PDBs\n";
    if (nPdbs == 0) return;

    // Dihedrals and side-chain burial come from NativeQCInputs.h.
    long long nChainBreaks = 0, nParamFallback = 0;
    cerr << "structure_survey: dihedrals and burial computed from coordinates\n";

    // ── Set up SurfacePatcher / CavityDetector ───────────────────────────
    SurfacePatcher patcher;
    patcher.burial_radius    = Option.surface_burial_radius;
    patcher.burial_threshold = Option.surface_burial_threshold;
    patcher.adj_cutoff       = Option.surface_adj_cutoff;
    patcher.min_patch_size   = Option.surface_min_patch_size;
    patcher.hphob_set        = Option.surface_hphob_set;
    patcher.pH               = Option.pH;

    CavityDetector cavDet;
    cavDet.grid_spacing = Option.cavity_grid_spacing;
    cavDet.probe_radius = Option.cavity_probe_radius;
    cavDet.min_volume   = Option.cavity_min_volume;

    // ── Reference DB ─────────────────────────────────────────────────────
    StructureReferenceDB refDB;

    // ── Open TSV output ──────────────────────────────────────────────────
    string tsvPath = Option.prefix + "structure_survey.tsv";
    ofstream tsv(tsvPath.c_str());
    // NB: `sc_burial` is FoldX-style side-chain burial (HIGH = buried), not an
    // accessibility. Renamed from the former `facc` header, which invited
    // downstream code to read it with the sign reversed. The energy columns
    // are kept for format compatibility and are left empty.
    tsv << "source_pdb\tchain\taa1\tresnum\tss\tsc_burial"
        << "\tomega\tphi\tpsi"
        << "\ttotal_energy\tbackbone_hbond\tsidechain_hbond"
        << "\tvdw\telectrostatics\tsolv_polar\tsolv_hydrophobic"
        << "\tvdw_clashes\tentropy_sc\tentropy_mc"
        << "\tsloop_entropy\tmloop_entropy"
        << "\tcis_bond\ttorsional_clash\tbackbone_clash"
        << "\thelix_dipole\twater_bridge\tdisulfide"
        << "\telectrostatic_kon\tpartial_covalent\tenergy_ionisation"
        << "\tentropy_complex\thetero_back_hbond\thetero_side_hbond"
        << "\tclashscore"
        << "\t" << BondGeometry::Header()
        << "\n";

    // ── Process in batches (bounds the number of PDBs held in memory) ─────
    int batchSize = 1000;
    auto t0 = chrono::steady_clock::now();

    // Open cavity output file
    string cavPath = Option.prefix + "structure_survey_cavities.tsv";
    ofstream cavOut(cavPath.c_str());
    cavOut << "pdb_name\tn_cavities\ttotal_cavity_vol\tlargest_cavity_vol"
           << "\tprotein_volume\tvoid_fraction\n";

    // Open patches output file
    string patchPath = Option.prefix + "structure_survey_patches.tsv";
    ofstream patchOut(patchPath.c_str());
    patchOut << "pdb_name\tchain\tn_hphob_patches\thphob1_size\thphob1_kd_sum"
             << "\tn_pos_clusters\tpos1_size\tpos1_charge"
             << "\tn_neg_clusters\tneg1_size\tneg1_charge\n";

    int nProcessed = 0, nFailed = 0;
    // Residues discarded because their residue name did not map to a
    // one-letter code. Must stay at zero: a non-zero count means rows are being
    // dropped from the per-residue table and the reference DB without notice.
    long long nUnmappedRes = 0;
    map<string, long long> unmappedNames;

    for (int batch0 = 0; batch0 < nPdbs; batch0 += batchSize)
    {
        int batchEnd = min(batch0 + batchSize, nPdbs);
        int batchN = batchEnd - batch0;

        // ── Read PDB files ───────────────────────────────────────────
        vector<PDB_File> pdbFiles(batchN);
        vector<string>   pdbNames(batchN);
        vector<bool>     pdbOk(batchN, false);

        for (int bi = 0; bi < batchN; bi++)
        {
            const string& path = pdbPaths[batch0 + bi];

            // Extract base name
            string baseName = path;
            size_t slash = baseName.find_last_of("/\\");
            if (slash != string::npos) baseName = baseName.substr(slash + 1);
            size_t dot = baseName.find_last_of('.');
            if (dot != string::npos) baseName = baseName.substr(0, dot);
            pdbNames[bi] = baseName;

            // Read PDB for our own analysis
            PDB_Reader reader;
            if (reader.Readpdb(path, pdbFiles[bi]) < 0) {
                cerr << "  SKIP " << baseName << " (read failed)\n";
                nFailed++;
                continue;
            }
            pdbOk[bi] = true;
        }

        // ── Process each PDB in the batch ────────────────────────────
        for (int bi = 0; bi < batchN; bi++)
        {
            if (!pdbOk[bi]) continue;
            const string& name = pdbNames[bi];
            const PDB_File& pdb = pdbFiles[bi];

            // Per-residue dihedrals + side-chain burial from the coordinates
            NativeQCSummary nqs;
            vector<NativeResidueInputs> nat = ComputeNativeQCInputs(pdb, &nqs);
            nChainBreaks   += nqs.nChainBreaks;
            nParamFallback += nqs.nParamFallback;
            if (nat.empty()) {
                cerr << "  SKIP " << name << " (no residues)\n";
                nFailed++;
                continue;
            }

            vector<SDParsed> sdData;
            sdData.reserve(nat.size());
            for (const auto& q : nat) {
                SDParsed sd;
                for (int m = 0; m < SM_COUNT; m++) sd.values[m] = NAN;
                sd.residueName = q.resName;
                sd.chain       = q.chain;
                sd.resNumber   = q.resnum;
                sd.ss          = '-';           // not assigned
                sd.facc        = q.scBurial;    // side-chain burial (HIGH = buried)
                sd.values[SM_OMEGA] = q.omega;
                sd.values[SM_PHI]   = q.phi;
                sd.values[SM_PSI]   = q.psi;
                sd.valid = true;
                sdData.push_back(sd);
            }

            // ── BondGeometry for each chain ──────────────────────────
            // Build a map: (chain, resnum) → BondGeometryResult
            unordered_map<long long, BondGeometryResult> bgMap;
            for (const auto& pep : pdb.peptides) {
                auto bgResults = BondGeometry::Compute(pep);
                for (const auto& bg : bgResults) {
                    long long key = ((long long)(unsigned char)bg.chain << 32) | (long long)bg.resnum;
                    bgMap[key] = bg;
                }
            }

            // ── MolProbity-style clash score (whole structure) ───────
            // Build a map: (chain, resnum) → per-residue clash count
            auto clashRes = Clashscore::Compute(pdb);
            unordered_map<long long, double> clashMap;
            for (const auto& cr : clashRes) {
                long long key = ((long long)(unsigned char)cr.chain << 32) | (long long)cr.resnum;
                clashMap[key] = cr.clashscore;
            }
            // heavy-atom count for the structure-level clashscore (per 1000 atoms)
            int nHeavyAtoms = 0;
            for (const auto& pep : pdb.peptides)
                for (const auto& res : pep.reslist) {
                    if (res.residue_name == "HOH" || res.residue_name == "GAP") continue;
                    for (const auto& a : res.atoms)
                        if (!a.atom_name.empty() && a.atom_name[0] != 'H') ++nHeavyAtoms;
                }

            // ── Merge SD + BondGeometry and accumulate into RefDB ────
            for (const auto& sd : sdData)
            {
                char aa1 = threeToOne(sd.residueName);
                if (aa1 == '?') {
                    nUnmappedRes++;
                    unmappedNames[sd.residueName]++;
                    continue;
                }

                // Look up bond geometry
                long long bgKey = ((long long)(unsigned char)sd.chain << 32) | (long long)sd.resNumber;
                double values[SM_COUNT];
                memcpy(values, sd.values, sizeof(values));

                // Clash score (whole-structure; default 0 where the residue had atoms)
                auto clIt = clashMap.find(bgKey);
                if (clIt != clashMap.end()) values[SM_CLASHSCORE] = clIt->second;

                auto bgIt = bgMap.find(bgKey);
                if (bgIt != bgMap.end()) {
                    const auto& bg = bgIt->second;
                    values[SM_BOND_N_CA]    = bg.bond_N_CA;
                    values[SM_BOND_CA_C]    = bg.bond_CA_C;
                    values[SM_BOND_C_N]     = bg.bond_C_N;
                    values[SM_BOND_C_O]     = bg.bond_C_O;
                    values[SM_BOND_CA_CB]   = bg.bond_CA_CB;
                    values[SM_ANGLE_N_CA_C] = bg.angle_N_CA_C;
                    values[SM_ANGLE_CA_C_N] = bg.angle_CA_C_N;
                    values[SM_ANGLE_C_N_CA] = bg.angle_C_N_CA;
                    values[SM_CB_DEVIATION] = bg.cb_deviation;
                    values[SM_CHI1]         = bg.chi1;
                    values[SM_CHI2]         = bg.chi2;
                    values[SM_CHI3]         = bg.chi3;
                    values[SM_CHI4]         = bg.chi4;
                    values[SM_PEPTIDE_PLANARITY] = bg.peptide_planarity;
                    values[SM_ROTAMER_OUTLIER]   = bg.rotamer_outlier;
                    values[SM_CABLAM_MU_IN]      = bg.cablam_mu_in;
                    values[SM_CABLAM_MU_OUT]     = bg.cablam_mu_out;
                    values[SM_CABLAM_NU]         = bg.cablam_nu;
                    values[SM_CABLAM_CA_VIRTUAL] = bg.cablam_ca_virtual;
                    values[SM_CIS_PEPTIDE]       = bg.cis_peptide;
                    values[SM_CA_CHIRALITY]      = bg.ca_chirality;
                }

                // Accumulate into reference DB
                refDB.AddObservation(aa1, sd.facc, values);

                // Accumulate Ramachandran (phi, psi) observation
                refDB.AddRamachandranObs(aa1, values[SM_PHI], values[SM_PSI]);

                // Write TSV row
                tsv << name << "\t" << sd.chain << "\t" << aa1
                    << "\t" << sd.resNumber << "\t" << sd.ss
                    << "\t" << fixed << setprecision(3) << sd.facc;
                for (int m = 0; m < SM_BOND_N_CA; m++)
                    tsv << "\t" << (isnan(values[m]) ? "" : to_string(values[m]));
                // Bond geometry columns
                if (bgIt != bgMap.end())
                    tsv << "\t" << BondGeometry::ToTSV(bgIt->second);
                else
                    tsv << "\t\t\t\t\t\t\t\t";
                tsv << "\n";
            }

            // ── Cavities ─────────────────────────────────────────────
            CavitySummary cav = cavDet.Detect(pdb);
            cavOut << name
                   << "\t" << cav.n_cavities
                   << "\t" << fixed << setprecision(1) << cav.total_volume
                   << "\t" << fixed << setprecision(1) << cav.largest_volume
                   << "\t" << fixed << setprecision(1) << cav.protein_volume
                   << "\t" << fixed << setprecision(4) << cav.void_fraction
                   << "\n";

            // ── Surface patches (per-chain, aggregated for whole-struct) ─
            int    ws_nHphob = 0, ws_nPos = 0, ws_nNeg = 0;
            int    ws_largestHphobSize = 0, ws_largestPosSize = 0, ws_largestNegSize = 0;
            double ws_largestHphobKd = 0, ws_largestPosCharge = 0, ws_largestNegCharge = 0;

            for (const auto& pep : pdb.peptides) {
                if (pep.reslist.empty()) continue;
                SurfaceAnalysis sa = patcher.Analyse(pep, name);
                patchOut << name << "\t" << sa.chain
                         << "\t" << sa.hydrophobic.size()
                         << "\t" << sa.topHydrophobic().size
                         << "\t" << fixed << setprecision(1) << sa.topHydrophobic().kd_sum
                         << "\t" << sa.positive.size()
                         << "\t" << sa.topPositive().size
                         << "\t" << fixed << setprecision(1) << sa.topPositive().net_charge
                         << "\t" << sa.negative.size()
                         << "\t" << sa.topNegative().size
                         << "\t" << fixed << setprecision(1) << sa.topNegative().net_charge
                         << "\n";

                // Aggregate across chains for whole-struct reference
                ws_nHphob += (int)sa.hydrophobic.size();
                ws_nPos   += (int)sa.positive.size();
                ws_nNeg   += (int)sa.negative.size();
                if (sa.topHydrophobic().size > ws_largestHphobSize) {
                    ws_largestHphobSize = sa.topHydrophobic().size;
                    ws_largestHphobKd   = sa.topHydrophobic().kd_sum;
                }
                if (sa.topPositive().size > ws_largestPosSize) {
                    ws_largestPosSize   = sa.topPositive().size;
                    ws_largestPosCharge = sa.topPositive().net_charge;
                }
                if (sa.topNegative().size > ws_largestNegSize) {
                    ws_largestNegSize   = sa.topNegative().size;
                    ws_largestNegCharge = fabs(sa.topNegative().net_charge);
                }
            }

            // ── Accumulate whole-structure observation ───────────────
            int nResTotal = (int)sdData.size();
            double wsValues[WS_COUNT];
            wsValues[WS_N_CAVITIES]         = cav.n_cavities;
            wsValues[WS_TOTAL_CAVITY_VOL]   = cav.total_volume;
            wsValues[WS_LARGEST_CAVITY_VOL] = cav.largest_volume;
            wsValues[WS_VOID_FRACTION]      = cav.void_fraction;
            wsValues[WS_N_HPHOB_PATCHES]    = ws_nHphob;
            wsValues[WS_LARGEST_HPHOB_SIZE] = ws_largestHphobSize;
            wsValues[WS_LARGEST_HPHOB_KD]   = ws_largestHphobKd;
            wsValues[WS_N_POS_CLUSTERS]     = ws_nPos;
            wsValues[WS_LARGEST_POS_SIZE]   = ws_largestPosSize;
            wsValues[WS_LARGEST_POS_CHARGE] = ws_largestPosCharge;
            wsValues[WS_N_NEG_CLUSTERS]     = ws_nNeg;
            wsValues[WS_LARGEST_NEG_SIZE]   = ws_largestNegSize;
            wsValues[WS_LARGEST_NEG_CHARGE] = ws_largestNegCharge;
            wsValues[WS_CLASHSCORE]         = Clashscore::StructureClashscore(clashRes, nHeavyAtoms);

            refDB.AddWholeStructObs(nResTotal, wsValues);

            nProcessed++;
            refDB.totalPDBs++;

            // ── Progress every 100 PDBs ─────────────────────────────
            if (nProcessed % 100 == 0 || nProcessed == nPdbs)
            {
                auto now = chrono::steady_clock::now();
                int elapsed = (int)chrono::duration_cast<chrono::seconds>(now - t0).count();
                double rate = nProcessed > 0 ? (double)elapsed / nProcessed : 0;
                int remaining = (int)(rate * (nPdbs - nProcessed));
                int etaH = remaining / 3600;
                int etaM = (remaining % 3600) / 60;
                cerr << "[survey] " << nProcessed << "/" << nPdbs
                     << " done (" << nFailed << " failed)"
                     << "  elapsed " << elapsed/60 << "m"
                     << "  ETA " << etaH << "h" << setfill('0') << setw(2) << etaM << "m"
                     << setfill(' ')
                     << "  (" << refDB.totalResidues << " residues)"
                     << endl;
            }
        }
    }

    cerr << "\n";

    // ── Finalise and write reference DB ──────────────────────────────────
    refDB.Finalise();

    string binPath = Option.prefix + "structure_reference.bin";
    if (refDB.WriteBinary(binPath))
        cerr << "structure_survey: reference DB written to " << binPath << "\n";
    else
        cerr << "structure_survey: ERROR writing " << binPath << "\n";

    // Human-readable summary
    string sumPath = Option.prefix + "structure_survey_summary.txt";
    ofstream sumOut(sumPath.c_str());
    refDB.PrintSummary(sumOut);

    cerr << "structure_survey: done. " << nProcessed << " PDBs, "
         << refDB.totalResidues << " residues.\n";
    if (nUnmappedRes > 0) {
        cerr << "  WARNING: " << nUnmappedRes << " residues dropped, "
             << "unrecognised residue name:";
        for (const auto& kv : unmappedNames)
            cerr << " " << kv.first << "(" << kv.second << ")";
        cerr << "\n           These are missing from the per-residue table "
             << "and the reference DB.\n";
    }
    if (nChainBreaks > 0)
        cerr << "  note: " << nChainBreaks << " chain breaks (C-N > 2.0 A); "
             << "the dihedrals spanning them were not added to the reference.\n";
    if (nParamFallback > 0)
        cerr << "  WARNING: " << nParamFallback << " atoms have no volume "
             << "parameter (generic values used for their burial).\n";
    cerr << "  reference DB: " << binPath << "\n";
    cerr << "  per-residue:  " << tsvPath << "\n";
    cerr << "  cavities:     " << cavPath << "\n";
    cerr << "  patches:      " << patchPath << "\n";
    cerr << "  summary:      " << sumPath << "\n";
}
