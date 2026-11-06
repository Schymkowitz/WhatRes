// pdb_quality_check_command.cpp
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz

#include "pdb_quality_check_command.h"
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
#include <algorithm>
#include <numeric>
#include <thread>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <map>

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
// Per-residue input record. Filled from NativeQCInputs.h (dihedrals and
// FoldX-style side-chain burial computed from the coordinates). The energy
// metrics stay NaN, so the scoring loop skips them.
// ─────────────────────────────────────────────────────────────────────────────

struct QCParsed
{
    string residueName;  // 3-letter
    char   chain;
    int    resNumber;
    char   ss;
    double facc;
    double values[SM_COUNT];
    bool   valid = false;
};

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

// ─────────────────────────────────────────────────────────────────────────────
// Scoring function
//
// For a single metric value and its reference distribution:
//   1. Find where the value falls in the percentile distribution
//   2. Convert to a "surprise" score:
//      score = -log2( min(pct, 100-pct) / 50 )
//
//   This gives: 0.0 at median, ~1.0 at p25/p75, ~3.3 at p5/p95,
//               ~5.6 at p1/p99, and higher for more extreme outliers.
//
// The intuition: doubling the rarity adds +1 to the score, so a value at
// the 1st percentile scores ~5.6 (about 50× rarer than median).
// ─────────────────────────────────────────────────────────────────────────────

static double metricScore(double value, const PercentileStats& ref)
{
    // If we have no distribution data, can't score
    if (ref.sd == 0.0 && ref.mean == 0.0) return 0.0;

    double pct = ref.Percentile(value);

    // Tail probability: how far from the middle?
    double tail = min(pct, 100.0 - pct);  // 0..50
    if (tail < 0.1) tail = 0.1;           // cap at p0.1 to avoid infinity

    // Surprise: 0 at median (tail=50), increases as tail shrinks
    double surprise = -log2(tail / 50.0);
    if (surprise < 0.0) surprise = 0.0;

    return surprise;
}

// Per-residue result
struct ResidueQC
{
    char   chain;
    char   aa1;
    int    resnum;
    char   ss;
    double facc;
    int    faccBin;
    int    refCount;              // observations in reference cell

    double metricValues[SM_COUNT];
    double metricPctile[SM_COUNT]; // percentile in reference distribution
    double metricScores[SM_COUNT]; // surprise scores
    bool   metricOutlier[SM_COUNT]; // outside p1–p99

    // Ramachandran 2D score (separate from 1D phi/psi scores)
    double ramaPctile;            // density-based percentile [0–100]
    double ramaScore;             // surprise score from 2D distribution
    bool   ramaOutlier;           // below threshold density

    double totalScore;            // sum of metric scores + rama score
    int    nOutliers;             // count of metrics outside p1–p99
    string worstMetric;           // name of highest-scoring metric
    double worstScore;
};

// Per-PDB summary
struct PDBQCSummary
{
    string pdbName;
    int    nResidues;
    int    nResiduesChecked;      // with valid reference cell
    double meanResidueScore;
    double maxResidueScore;
    double pdbScore;              // combined score
    int    totalOutliers;         // sum of per-residue outlier counts
    double worst5_mean;           // mean score of 5 worst residues

    // Whole-structure metrics
    double wsValues[WS_COUNT];
    double wsScores[WS_COUNT];
    double wsTotalScore;          // sum of whole-struct metric scores
    int    wsOutliers;            // whole-struct outlier count
    string wsWorstMetric;
    double wsWorstScore;

    // Final combined score
    double combinedScore;         // pdbScore + wsTotalScore contribution
};


// ═════════════════════════════════════════════════════════════════════════════
// Execute
// ═════════════════════════════════════════════════════════════════════════════
void pdb_quality_check_command::Execute()
{
    if (Option.inputfiles.empty())
    {
        cerr << "pdb_quality_check: no input file specified.\n";
        return;
    }

    // ── Load reference database ─────────────────────────────────────────
    string refPath = Option.qc_reference_db;
    if (refPath.empty()) {
        cerr << "pdb_quality_check: no reference DB specified.\n"
             << "  Set <qc_reference_db>structure_reference.bin</qc_reference_db> "
             << "in your Options file.\n";
        return;
    }

    StructureReferenceDB refDB;
    if (!refDB.ReadBinary(refPath)) {
        cerr << "pdb_quality_check: failed to read " << refPath << "\n";
        return;
    }
    cerr << "pdb_quality_check: loaded reference DB (" << refDB.totalPDBs
         << " PDBs, " << refDB.totalResidues << " residues)\n";

    // ── Scoring thresholds from options ──────────────────────────────────
    // Outlier percentile indices: default p1 (index 0) and p99 (index 6)
    int loIdx = Option.qc_outlier_lo_pctile;  // 0=p1, 1=p5, 2=p25
    int hiIdx = Option.qc_outlier_hi_pctile;  // 4=p75, 5=p95, 6=p99
    if (loIdx < 0 || loIdx > 6) loIdx = 0;
    if (hiIdx < 0 || hiIdx > 6) hiIdx = 6;

    // Minimum reference count to trust a cell
    int minRefCount = Option.qc_min_ref_count;
    if (minRefCount <= 0) minRefCount = 50;

    // ── Read PDB list ───────────────────────────────────────────────────
    vector<string> pdbPaths;
    {
        ifstream list(Option.inputfiles[0].c_str());
        if (!list) {
            cerr << "pdb_quality_check: cannot open " << Option.inputfiles[0] << endl;
            return;
        }
        string line;
        while (getline(list, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (!Option.surface_pdb_dir.empty() && line[0] != '/' && line[0] != '~') {
                string dir = Option.surface_pdb_dir;
                if (dir.back() != '/') dir += '/';
                line = dir + line;
            }
            pdbPaths.push_back(line);
        }
    }

    int nPdbs = (int)pdbPaths.size();
    cerr << "pdb_quality_check: " << nPdbs << " PDBs to check\n";
    if (nPdbs == 0) return;

    // Residues discarded because their residue name did not map to a
    // one-letter code. Must stay at zero: a non-zero count means residues are
    // being skipped from the QC without notice.
    long long nUnmappedRes = 0;
    map<string, long long> unmappedNames;

    // Dihedrals and side-chain burial come from NativeQCInputs.h.
    // Totals over all structures, reported at the end so nothing is dropped silently.
    long long nChainBreaks = 0, nParamFallback = 0;
    cerr << "pdb_quality_check: dihedrals and burial computed from coordinates\n";

    // ── Set up CavityDetector / SurfacePatcher ──────────────────────────
    CavityDetector cavDet;
    cavDet.grid_spacing = Option.cavity_grid_spacing;
    cavDet.probe_radius = Option.cavity_probe_radius;
    cavDet.min_volume   = Option.cavity_min_volume;

    SurfacePatcher patcher;
    patcher.burial_radius    = Option.surface_burial_radius;
    patcher.burial_threshold = Option.surface_burial_threshold;
    patcher.adj_cutoff       = Option.surface_adj_cutoff;
    patcher.min_patch_size   = Option.surface_min_patch_size;
    patcher.hphob_set        = Option.surface_hphob_set;
    patcher.pH               = Option.pH;

    // ── Open output files ───────────────────────────────────────────────
    string reportPath   = Option.prefix + "quality_report.tsv";
    string outlierPath  = Option.prefix + "quality_outliers.tsv";
    string summaryPath  = Option.prefix + "quality_summary.txt";

    ofstream report(reportPath.c_str());
    ofstream outliers(outlierPath.c_str());

    // Report header
    report << "pdb\tchain\taa\tresnum\tss\tsc_burial\tburial_bin\tref_count"
           << "\tresidue_score\tn_outliers\tworst_metric\tworst_score"
           << "\trama2D_pctile\trama2D_score";
    for (int m = 0; m < SM_COUNT; m++)
        report << "\t" << StructMetricName(m) << "_pctile"
               << "\t" << StructMetricName(m) << "_score";
    report << "\n";

    // Outlier header
    outliers << "pdb\tchain\taa\tresnum\tss\tsc_burial\tburial_bin"
             << "\tmetric\tvalue\tref_p1\tref_p50\tref_p99"
             << "\tpercentile\tscore\n";

    vector<PDBQCSummary> pdbSummaries;

    // ── Process PDBs ────────────────────────────────────────────────────
    int batchSize = 50;   // PDBs held in memory at once

    for (int batch0 = 0; batch0 < nPdbs; batch0 += batchSize)
    {
        int batchEnd = min(batch0 + batchSize, nPdbs);
        int batchN = batchEnd - batch0;

        // Read PDB files
        vector<PDB_File> pdbFiles(batchN);
        vector<string>   pdbNames(batchN);
        vector<bool>     pdbOk(batchN, false);

        for (int bi = 0; bi < batchN; bi++)
        {
            const string& path = pdbPaths[batch0 + bi];

            string baseName = path;
            size_t slash = baseName.find_last_of("/\\");
            if (slash != string::npos) baseName = baseName.substr(slash + 1);
            size_t dot = baseName.find_last_of('.');
            if (dot != string::npos) baseName = baseName.substr(0, dot);
            pdbNames[bi] = baseName;

            PDB_Reader reader;
            if (reader.Readpdb(path, pdbFiles[bi]) < 0) {
                cerr << "  SKIP " << baseName << " (read failed)\n";
                continue;
            }
            pdbOk[bi] = true;
        }

        // Process each PDB
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
            if (nat.empty()) { cerr << "  SKIP " << name << " (no residues)\n"; continue; }

            vector<QCParsed> sdData;
            sdData.reserve(nat.size());
            for (const auto& q : nat) {
                QCParsed sd;
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

            // BondGeometry
            unordered_map<long long, BondGeometryResult> bgMap;
            for (const auto& pep : pdb.peptides) {
                auto bgResults = BondGeometry::Compute(pep);
                for (const auto& bg : bgResults) {
                    long long key = ((long long)(unsigned char)bg.chain << 32) | (long long)bg.resnum;
                    bgMap[key] = bg;
                }
            }

            // MolProbity-style clash score (whole structure)
            auto clashRes = Clashscore::Compute(pdb);
            unordered_map<long long, double> clashMap;
            for (const auto& cr : clashRes) {
                long long key = ((long long)(unsigned char)cr.chain << 32) | (long long)cr.resnum;
                clashMap[key] = cr.clashscore;
            }
            int nHeavyAtoms = 0;
            for (const auto& pep : pdb.peptides)
                for (const auto& res : pep.reslist) {
                    if (res.residue_name == "HOH" || res.residue_name == "GAP") continue;
                    for (const auto& a : res.atoms)
                        if (!a.atom_name.empty() && a.atom_name[0] != 'H') ++nHeavyAtoms;
                }

            // ── Score each residue ──────────────────────────────────────
            vector<ResidueQC> residueResults;
            residueResults.reserve(sdData.size());

            for (const auto& sd : sdData)
            {
                char aa1 = threeToOne(sd.residueName);
                if (aa1 == '?') {
                    nUnmappedRes++;
                    unmappedNames[sd.residueName]++;
                    continue;
                }

                ResidueQC rqc;
                rqc.chain   = sd.chain;
                rqc.aa1     = aa1;
                rqc.resnum  = sd.resNumber;
                rqc.ss      = sd.ss;
                rqc.facc    = sd.facc;
                rqc.faccBin = FaccBin(sd.facc);

                // Merge SD + BondGeometry
                memcpy(rqc.metricValues, sd.values, sizeof(sd.values));
                long long bgKey = ((long long)(unsigned char)sd.chain << 32) | (long long)sd.resNumber;
                auto clIt = clashMap.find(bgKey);
                if (clIt != clashMap.end()) rqc.metricValues[SM_CLASHSCORE] = clIt->second;
                auto bgIt = bgMap.find(bgKey);
                if (bgIt != bgMap.end()) {
                    const auto& bg = bgIt->second;
                    rqc.metricValues[SM_BOND_N_CA]    = bg.bond_N_CA;
                    rqc.metricValues[SM_BOND_CA_C]    = bg.bond_CA_C;
                    rqc.metricValues[SM_BOND_C_N]     = bg.bond_C_N;
                    rqc.metricValues[SM_BOND_C_O]     = bg.bond_C_O;
                    rqc.metricValues[SM_BOND_CA_CB]   = bg.bond_CA_CB;
                    rqc.metricValues[SM_ANGLE_N_CA_C] = bg.angle_N_CA_C;
                    rqc.metricValues[SM_ANGLE_CA_C_N] = bg.angle_CA_C_N;
                    rqc.metricValues[SM_ANGLE_C_N_CA] = bg.angle_C_N_CA;
                    rqc.metricValues[SM_CB_DEVIATION] = bg.cb_deviation;
                    rqc.metricValues[SM_CHI1]         = bg.chi1;
                    rqc.metricValues[SM_CHI2]         = bg.chi2;
                    rqc.metricValues[SM_CHI3]         = bg.chi3;
                    rqc.metricValues[SM_CHI4]         = bg.chi4;
                    rqc.metricValues[SM_PEPTIDE_PLANARITY] = bg.peptide_planarity;
                    rqc.metricValues[SM_ROTAMER_OUTLIER]   = bg.rotamer_outlier;
                    rqc.metricValues[SM_CABLAM_MU_IN]      = bg.cablam_mu_in;
                    rqc.metricValues[SM_CABLAM_MU_OUT]     = bg.cablam_mu_out;
                    rqc.metricValues[SM_CABLAM_NU]         = bg.cablam_nu;
                    rqc.metricValues[SM_CABLAM_CA_VIRTUAL] = bg.cablam_ca_virtual;
                    rqc.metricValues[SM_CIS_PEPTIDE]       = bg.cis_peptide;
                    rqc.metricValues[SM_CA_CHIRALITY]      = bg.ca_chirality;
                }

                // Look up reference cell
                const RefCell* cell = refDB.GetCell(aa1, sd.facc);
                rqc.refCount = cell ? cell->count : 0;

                // Ramachandran 2D scoring
                rqc.ramaPctile  = -1.0;
                rqc.ramaScore   = 0.0;
                rqc.ramaOutlier = false;

                if (refDB.HasRamachandranData() &&
                    !isnan(rqc.metricValues[SM_PHI]) && !isnan(rqc.metricValues[SM_PSI]))
                {
                    double rp = refDB.GetRamachandranPctile(
                        aa1, rqc.metricValues[SM_PHI], rqc.metricValues[SM_PSI]);
                    if (rp >= 0.0)
                    {
                        rqc.ramaPctile = rp;
                        // Treat density percentile like a 1D percentile for scoring:
                        // high percentile = common region = low surprise
                        // We score based on how rare the region is:
                        //   score = -log2( max(pctile, 0.2) / 50 )
                        double tail = std::max(rp, 0.2);  // cap to avoid infinity
                        rqc.ramaScore = -log2(tail / 50.0);
                        if (rqc.ramaScore < 0.0) rqc.ramaScore = 0.0;
                        // Outlier if in bottom 2% of density (pctile < 2)
                        rqc.ramaOutlier = (rp < 2.0);
                    }
                }

                // Score each metric
                rqc.totalScore  = rqc.ramaScore;  // start with rama contribution
                rqc.nOutliers   = rqc.ramaOutlier ? 1 : 0;
                rqc.worstScore  = rqc.ramaScore;
                rqc.worstMetric = rqc.ramaScore > 0.0 ? "rama2D" : "";

                for (int m = 0; m < SM_COUNT; m++)
                {
                    rqc.metricPctile[m] = 50.0;
                    rqc.metricScores[m] = 0.0;
                    rqc.metricOutlier[m] = false;

                    if (isnan(rqc.metricValues[m])) continue;
                    if (!cell || cell->count < minRefCount) continue;

                    const PercentileStats& ref = cell->stats[m];
                    rqc.metricPctile[m] = ref.Percentile(rqc.metricValues[m]);
                    rqc.metricScores[m] = metricScore(rqc.metricValues[m], ref);
                    rqc.metricOutlier[m] = ref.IsOutlier(rqc.metricValues[m], loIdx, hiIdx);

                    rqc.totalScore += rqc.metricScores[m];
                    if (rqc.metricOutlier[m]) rqc.nOutliers++;

                    if (rqc.metricScores[m] > rqc.worstScore) {
                        rqc.worstScore = rqc.metricScores[m];
                        rqc.worstMetric = StructMetricName(m);
                    }

                    // Write outlier detail line
                    if (rqc.metricOutlier[m]) {
                        outliers << name << "\t" << sd.chain << "\t" << aa1
                                 << "\t" << sd.resNumber << "\t" << sd.ss
                                 << "\t" << fixed << setprecision(3) << sd.facc
                                 << "\t" << FaccBinLabel(rqc.faccBin)
                                 << "\t" << StructMetricName(m)
                                 << "\t" << fixed << setprecision(3) << rqc.metricValues[m]
                                 << "\t" << fixed << setprecision(3) << ref.p[0]    // p1
                                 << "\t" << fixed << setprecision(3) << ref.p[3]    // p50
                                 << "\t" << fixed << setprecision(3) << ref.p[6]    // p99
                                 << "\t" << fixed << setprecision(1) << rqc.metricPctile[m]
                                 << "\t" << fixed << setprecision(2) << rqc.metricScores[m]
                                 << "\n";
                    }
                }

                // Write Ramachandran outlier line
                if (rqc.ramaOutlier)
                {
                    outliers << name << "\t" << sd.chain << "\t" << aa1
                             << "\t" << sd.resNumber << "\t" << sd.ss
                             << "\t" << fixed << setprecision(3) << sd.facc
                             << "\t" << FaccBinLabel(rqc.faccBin)
                             << "\trama2D"
                             << "\t" << fixed << setprecision(1) << rqc.ramaPctile
                             << "\t-\t-\t-"
                             << "\t" << fixed << setprecision(1) << rqc.ramaPctile
                             << "\t" << fixed << setprecision(2) << rqc.ramaScore
                             << "\n";
                }

                // Write full report line
                report << name << "\t" << rqc.chain << "\t" << rqc.aa1
                       << "\t" << rqc.resnum << "\t" << rqc.ss
                       << "\t" << fixed << setprecision(3) << rqc.facc
                       << "\t" << FaccBinLabel(rqc.faccBin)
                       << "\t" << rqc.refCount
                       << "\t" << fixed << setprecision(2) << rqc.totalScore
                       << "\t" << rqc.nOutliers
                       << "\t" << rqc.worstMetric
                       << "\t" << fixed << setprecision(2) << rqc.worstScore
                       << "\t" << fixed << setprecision(1) << rqc.ramaPctile
                       << "\t" << fixed << setprecision(2) << rqc.ramaScore;

                for (int m = 0; m < SM_COUNT; m++)
                    report << "\t" << fixed << setprecision(1) << rqc.metricPctile[m]
                           << "\t" << fixed << setprecision(2) << rqc.metricScores[m];
                report << "\n";

                residueResults.push_back(rqc);
            }

            // ── PDB-level summary ───────────────────────────────────────
            PDBQCSummary psum;
            psum.pdbName          = name;
            psum.nResidues        = (int)sdData.size();
            psum.nResiduesChecked = (int)residueResults.size();
            psum.totalOutliers    = 0;

            // Compute mean and max residue scores
            double sumScore = 0.0;
            psum.maxResidueScore = 0.0;
            vector<double> allScores;
            allScores.reserve(residueResults.size());

            for (const auto& rqc : residueResults) {
                sumScore += rqc.totalScore;
                allScores.push_back(rqc.totalScore);
                psum.totalOutliers += rqc.nOutliers;
                if (rqc.totalScore > psum.maxResidueScore)
                    psum.maxResidueScore = rqc.totalScore;
            }

            psum.meanResidueScore = allScores.empty() ? 0.0 : sumScore / allScores.size();

            // Worst-5 mean: average of top 5 residue scores
            sort(allScores.begin(), allScores.end(), greater<double>());
            int nWorst = min(5, (int)allScores.size());
            double w5sum = 0.0;
            for (int i = 0; i < nWorst; i++) w5sum += allScores[i];
            psum.worst5_mean = nWorst > 0 ? w5sum / nWorst : 0.0;

            // Per-residue PDB score
            psum.pdbScore = psum.meanResidueScore + 0.5 * psum.worst5_mean;

            // ── Whole-structure scoring (cavities + patches) ────────
            CavitySummary cav = cavDet.Detect(pdb);

            // Aggregate patches across chains (same logic as survey)
            int    ws_nHphob = 0, ws_nPos = 0, ws_nNeg = 0;
            int    ws_largestHphobSize = 0, ws_largestPosSize = 0, ws_largestNegSize = 0;
            double ws_largestHphobKd = 0, ws_largestPosCharge = 0, ws_largestNegCharge = 0;

            for (const auto& pep : pdb.peptides) {
                if (pep.reslist.empty()) continue;
                SurfaceAnalysis sa = patcher.Analyse(pep, name);
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

            // Fill whole-structure values
            psum.wsValues[WS_N_CAVITIES]         = cav.n_cavities;
            psum.wsValues[WS_TOTAL_CAVITY_VOL]   = cav.total_volume;
            psum.wsValues[WS_LARGEST_CAVITY_VOL] = cav.largest_volume;
            psum.wsValues[WS_VOID_FRACTION]      = cav.void_fraction;
            psum.wsValues[WS_N_HPHOB_PATCHES]    = ws_nHphob;
            psum.wsValues[WS_LARGEST_HPHOB_SIZE] = ws_largestHphobSize;
            psum.wsValues[WS_LARGEST_HPHOB_KD]   = ws_largestHphobKd;
            psum.wsValues[WS_N_POS_CLUSTERS]     = ws_nPos;
            psum.wsValues[WS_LARGEST_POS_SIZE]   = ws_largestPosSize;
            psum.wsValues[WS_LARGEST_POS_CHARGE] = ws_largestPosCharge;
            psum.wsValues[WS_N_NEG_CLUSTERS]     = ws_nNeg;
            psum.wsValues[WS_LARGEST_NEG_SIZE]   = ws_largestNegSize;
            psum.wsValues[WS_LARGEST_NEG_CHARGE] = ws_largestNegCharge;
            psum.wsValues[WS_CLASHSCORE]         = Clashscore::StructureClashscore(clashRes, nHeavyAtoms);

            // Score against reference whole-structure distribution
            const WholeStructCell* wsCell = refDB.GetWSCell(psum.nResidues);
            psum.wsTotalScore = 0.0;
            psum.wsOutliers   = 0;
            psum.wsWorstScore = 0.0;
            psum.wsWorstMetric = "";

            for (int m = 0; m < WS_COUNT; m++)
            {
                psum.wsScores[m] = 0.0;
                if (!wsCell || wsCell->count < minRefCount) continue;

                const PercentileStats& ref = wsCell->stats[m];
                psum.wsScores[m] = metricScore(psum.wsValues[m], ref);
                psum.wsTotalScore += psum.wsScores[m];

                if (ref.IsOutlier(psum.wsValues[m], loIdx, hiIdx)) {
                    psum.wsOutliers++;
                    // Write to outliers file
                    outliers << name << "\t-\t-\t-\t-\t-\t"
                             << SizeBinLabel(SizeBin(psum.nResidues))
                             << "\tWS:" << WholeStructMetricName(m)
                             << "\t" << fixed << setprecision(3) << psum.wsValues[m]
                             << "\t" << fixed << setprecision(3) << ref.p[0]
                             << "\t" << fixed << setprecision(3) << ref.p[3]
                             << "\t" << fixed << setprecision(3) << ref.p[6]
                             << "\t" << fixed << setprecision(1) << ref.Percentile(psum.wsValues[m])
                             << "\t" << fixed << setprecision(2) << psum.wsScores[m]
                             << "\n";
                }

                if (psum.wsScores[m] > psum.wsWorstScore) {
                    psum.wsWorstScore = psum.wsScores[m];
                    psum.wsWorstMetric = WholeStructMetricName(m);
                }
            }

            // Combined score: per-residue + whole-structure contribution
            // Weight whole-struct at 0.25 to keep per-residue dominant
            psum.combinedScore = psum.pdbScore + 0.25 * psum.wsTotalScore;

            pdbSummaries.push_back(psum);

            cerr << "  " << name
                 << ": combined=" << fixed << setprecision(2) << psum.combinedScore
                 << "  res_score=" << psum.pdbScore
                 << "  ws_score=" << psum.wsTotalScore
                 << "  outliers=" << psum.totalOutliers
                 << "+" << psum.wsOutliers << "ws"
                 << "  residues=" << psum.nResiduesChecked << "\n";
        }
    }

    if (nUnmappedRes > 0) {
        cerr << "  WARNING: " << nUnmappedRes << " residues skipped, "
             << "unrecognised residue name:";
        for (const auto& kv : unmappedNames)
            cerr << " " << kv.first << "(" << kv.second << ")";
        cerr << "\n           These residues received no quality score.\n";
    }
    if (nChainBreaks > 0)
        cerr << "  note: " << nChainBreaks << " chain breaks (C-N > 2.0 A); "
             << "the dihedrals spanning them were left unscored.\n";
    if (nParamFallback > 0)
        cerr << "  WARNING: " << nParamFallback << " atoms have no volume "
             << "parameter (generic values used for their burial).\n";

    // ── Write summary file ──────────────────────────────────────────────
    ofstream summary(summaryPath.c_str());
    summary << "PDB Quality Check Summary\n"
            << "Reference DB: " << refPath << " ("
            << refDB.totalPDBs << " PDBs, " << refDB.totalResidues << " residues)\n"
            << "Outlier thresholds: p" << (loIdx == 0 ? 1 : loIdx == 1 ? 5 : 25)
            << " – p" << (hiIdx == 6 ? 99 : hiIdx == 5 ? 95 : 75) << "\n"
            << "Min reference count: " << minRefCount << "\n\n";

    // Sort by combined score (worst first)
    sort(pdbSummaries.begin(), pdbSummaries.end(),
         [](const PDBQCSummary& a, const PDBQCSummary& b) {
             return a.combinedScore > b.combinedScore;
         });

    summary << "PDB\tcombined_score\tres_score\tws_score\tmean_res_score\tworst5_mean"
            << "\tn_residues\tres_outliers\tws_outliers\tws_worst_metric";
    for (int m = 0; m < WS_COUNT; m++)
        summary << "\t" << WholeStructMetricName(m);
    summary << "\n";

    for (const auto& ps : pdbSummaries)
    {
        summary << ps.pdbName
                << "\t" << fixed << setprecision(2) << ps.combinedScore
                << "\t" << fixed << setprecision(2) << ps.pdbScore
                << "\t" << fixed << setprecision(2) << ps.wsTotalScore
                << "\t" << fixed << setprecision(2) << ps.meanResidueScore
                << "\t" << fixed << setprecision(2) << ps.worst5_mean
                << "\t" << ps.nResiduesChecked
                << "\t" << ps.totalOutliers
                << "\t" << ps.wsOutliers
                << "\t" << ps.wsWorstMetric;
        for (int m = 0; m < WS_COUNT; m++)
            summary << "\t" << fixed << setprecision(2) << ps.wsValues[m];
        summary << "\n";
    }

    // Overall statistics
    if (!pdbSummaries.empty()) {
        double totalCombined = 0.0;
        for (const auto& ps : pdbSummaries) totalCombined += ps.combinedScore;
        double mean = totalCombined / pdbSummaries.size();

        summary << "\n--- Overall ---\n"
                << "PDBs checked: " << pdbSummaries.size() << "\n"
                << "Mean combined score: " << fixed << setprecision(2) << mean << "\n"
                << "Worst PDB: " << pdbSummaries.front().pdbName
                << " (combined=" << fixed << setprecision(2) << pdbSummaries.front().combinedScore << ")\n"
                << "Best PDB: " << pdbSummaries.back().pdbName
                << " (combined=" << fixed << setprecision(2) << pdbSummaries.back().combinedScore << ")\n";
    }

    cerr << "\npdb_quality_check: done.\n"
         << "  report:   " << reportPath << "\n"
         << "  outliers: " << outlierPath << "\n"
         << "  summary:  " << summaryPath << "\n";
}
