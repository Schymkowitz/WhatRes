// StructureReferenceDB.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Reference database of structural metric distributions, built from a survey
// of ~25k high-resolution PDB structures.
//
// The database is a 2D lookup table:  amino acid type (20) × Facc bin (5).
// Each cell stores percentile arrays for every metric (dihedrals, energies,
// bond lengths, bond angles), plus rotamer population counts.
//
// The survey command accumulates raw values, then calls Finalise() to compute
// percentiles.  The result is serialised to a compact binary file.
//
// The quality check command loads the binary, then for each residue in an
// input PDB, looks up the appropriate cell and checks whether the residue's
// metrics fall within the expected range.
//
// Binary format:
//   [magic: "SRDB"]  [version: uint32]
//   [nAA: uint32]    [nFaccBins: uint32]  [nMetrics: uint32]
//   For each AA × Facc cell:
//     [count: uint32]
//     For each metric:
//       [7 percentiles: double × 7]   (p1, p5, p25, p50, p75, p95, p99)
//       [mean: double]  [sd: double]
//   [whole-structure distributions: cavity, patches, etc. (future extension)]
//
// Usage (survey):
//   StructureReferenceDB db;
//   db.AddObservation('V', 0.05, metric_values);  // buried Val
//   db.AddObservation('K', 0.65, metric_values);  // exposed Lys
//   ...
//   db.Finalise();
//   db.WriteBinary("structure_reference.bin");
//
// Usage (checker):
//   StructureReferenceDB db;
//   db.ReadBinary("structure_reference.bin");
//   auto* cell = db.GetCell('V', 0.05);
//   bool outlier = cell->IsOutlier(METRIC_PHI, measured_phi, 0.01, 0.99);
//
// -----------------------------------------------------------------------------

#pragma once
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <mutex>

// ─────────────────────────────────────────────────────────────────────────────
// Metric indices: these match the order in the per-residue survey output.
// FoldX SequenceDetail columns + BondGeometry columns.
// ─────────────────────────────────────────────────────────────────────────────
enum StructMetric
{
    // FoldX SequenceDetail dihedrals
    SM_OMEGA = 0, SM_PHI, SM_PSI,
    // FoldX SequenceDetail energies
    SM_TOTAL_ENERGY, SM_BACKBONE_HBOND, SM_SIDECHAIN_HBOND,
    SM_VDW, SM_ELECTROSTATICS, SM_SOLV_POLAR, SM_SOLV_HYDROPHOBIC,
    SM_VDW_CLASHES, SM_ENTROPY_SC, SM_ENTROPY_MC,
    SM_SLOOP_ENTROPY, SM_MLOOP_ENTROPY,
    SM_CIS_BOND, SM_TORSIONAL_CLASH, SM_BACKBONE_CLASH,
    SM_HELIX_DIPOLE, SM_WATER_BRIDGE, SM_DISULFIDE,
    SM_ELECTROSTATIC_KON, SM_PARTIAL_COVALENT, SM_ENERGY_IONISATION,
    SM_ENTROPY_COMPLEX,
    SM_HETERO_BACK_HBOND, SM_HETERO_SIDE_HBOND,
    // Coordinate-geometry validator: MolProbity-style all-atom clash count
    SM_CLASHSCORE,
    // BondGeometry
    SM_BOND_N_CA, SM_BOND_CA_C, SM_BOND_C_N, SM_BOND_C_O, SM_BOND_CA_CB,
    SM_ANGLE_N_CA_C, SM_ANGLE_CA_C_N, SM_ANGLE_C_N_CA,
    // Coordinate-geometry validators (from BondGeometry)
    SM_CB_DEVIATION,
    // Sidechain dihedrals (per amino-acid type)
    SM_CHI1, SM_CHI2, SM_CHI3, SM_CHI4,
    // Peptide-bond planarity and rotamer-outlier proxy
    SM_PEPTIDE_PLANARITY, SM_ROTAMER_OUTLIER,
    // CaBLAM backbone measures, cis-peptide, Cα chirality
    SM_CABLAM_MU_IN, SM_CABLAM_MU_OUT, SM_CABLAM_NU, SM_CABLAM_CA_VIRTUAL,
    SM_CIS_PEPTIDE, SM_CA_CHIRALITY,

    SM_COUNT   // total number of metrics
};

[[maybe_unused]] static const char* StructMetricName(int m)
{
    static const char* names[] = {
        "omega", "phi", "psi",
        "total_energy", "backbone_hbond", "sidechain_hbond",
        "vdw", "electrostatics", "solv_polar", "solv_hydrophobic",
        "vdw_clashes", "entropy_sc", "entropy_mc",
        "sloop_entropy", "mloop_entropy",
        "cis_bond", "torsional_clash", "backbone_clash",
        "helix_dipole", "water_bridge", "disulfide",
        "electrostatic_kon", "partial_covalent", "energy_ionisation",
        "entropy_complex",
        "hetero_back_hbond", "hetero_side_hbond",
        "clashscore",
        "bond_N_CA", "bond_CA_C", "bond_C_N", "bond_C_O", "bond_CA_CB",
        "angle_N_CA_C", "angle_CA_C_N", "angle_C_N_CA",
        "cb_deviation",
        "chi1", "chi2", "chi3", "chi4",
        "peptide_planarity", "rotamer_outlier",
        "cablam_mu_in", "cablam_mu_out", "cablam_nu", "cablam_ca_virtual",
        "cis_peptide", "ca_chirality"
    };
    return (m >= 0 && m < SM_COUNT) ? names[m] : "unknown";
}

// ─────────────────────────────────────────────────────────────────────────────
// Side-chain burial binning
//
// The binned quantity is FoldX's side-chain BURIAL (SequenceDetail col 32):
// HIGH = buried, LOW = exposed. It is not an accessibility. Note that glycine
// has no side chain and is therefore always exactly 0, so every glycine lands
// in bin 0 regardless of its actual environment.
// ─────────────────────────────────────────────────────────────────────────────
static const int    N_FACC_BINS = 5;
static const double FACC_EDGES[N_FACC_BINS + 1] = {0.0, 0.05, 0.20, 0.50, 0.80, 1.01};
//   bin 0: [0.00, 0.05)  fully exposed side chain (and all glycines)
//   bin 1: [0.05, 0.20)  mostly exposed
//   bin 2: [0.20, 0.50)  intermediate
//   bin 3: [0.50, 0.80)  mostly buried
//   bin 4: [0.80, 1.00]  fully buried

static int FaccBin(double burial)
{
    for (int b = 0; b < N_FACC_BINS; b++)
        if (burial < FACC_EDGES[b + 1]) return b;
    return N_FACC_BINS - 1;
}

static const char* FaccBinLabel(int b)
{
    static const char* labels[] = {
        "exposed_sc(0-5%)", "mostly_exposed(5-20%)", "intermediate(20-50%)",
        "mostly_buried(50-80%)", "buried(80-100%)"
    };
    return (b >= 0 && b < N_FACC_BINS) ? labels[b] : "?";
}

// ─────────────────────────────────────────────────────────────────────────────
// AA indexing (standard 20 amino acids)
// ─────────────────────────────────────────────────────────────────────────────
static const char AA_ORDER[] = "ACDEFGHIKLMNPQRSTVWY";
static const int  N_AA = 20;

static int AAIndex(char aa)
{
    for (int i = 0; i < N_AA; i++)
        if (AA_ORDER[i] == aa) return i;
    return -1;  // non-standard
}

// ─────────────────────────────────────────────────────────────────────────────
// PercentileStats: stats for one metric in one cell
// ─────────────────────────────────────────────────────────────────────────────
struct PercentileStats
{
    double p[7]  = {};     // p1, p5, p25, p50, p75, p95, p99
    double mean  = 0.0;
    double sd    = 0.0;

    // Check if a value is outside the [pLo, pHi] percentile range.
    // pLoIdx/pHiIdx: index into p[] array (0=p1, 1=p5, ..., 6=p99)
    bool IsOutlier(double value, int pLoIdx = 0, int pHiIdx = 6) const
    {
        return value < p[pLoIdx] || value > p[pHiIdx];
    }

    // Return which percentile this value falls at (approximate)
    double Percentile(double value) const
    {
        static const double pctValues[] = {1, 5, 25, 50, 75, 95, 99};
        if (value <= p[0]) return 0.5;    // below p1
        if (value >= p[6]) return 99.5;   // above p99
        for (int i = 0; i < 6; i++)
        {
            if (value >= p[i] && value <= p[i + 1])
            {
                if (std::abs(p[i + 1] - p[i]) < 1e-12) return pctValues[i];
                double frac = (value - p[i]) / (p[i + 1] - p[i]);
                return pctValues[i] + frac * (pctValues[i + 1] - pctValues[i]);
            }
        }
        return 50.0;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// RefCell: one AA × Facc bin cell
// ─────────────────────────────────────────────────────────────────────────────
struct RefCell
{
    int count = 0;                        // number of observations
    PercentileStats stats[SM_COUNT];      // one per metric

    // --- Accumulation phase (survey) ---
    std::vector<double> rawValues[SM_COUNT];  // raw values for percentile computation

    void AddObservation(const double values[SM_COUNT])
    {
        count++;
        for (int m = 0; m < SM_COUNT; m++)
            if (!std::isnan(values[m]))
                rawValues[m].push_back(values[m]);
    }

    void Finalise()
    {
        static const double pctTargets[] = {0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};

        for (int m = 0; m < SM_COUNT; m++)
        {
            auto& v = rawValues[m];
            int n = (int)v.size();
            if (n == 0) continue;

            std::sort(v.begin(), v.end());

            // Percentiles
            for (int pi = 0; pi < 7; pi++)
            {
                double idx = pctTargets[pi] * (n - 1);
                int lo = (int)std::floor(idx);
                int hi = std::min(lo + 1, n - 1);
                double frac = idx - lo;
                stats[m].p[pi] = v[lo] * (1.0 - frac) + v[hi] * frac;
            }

            // Mean & SD
            double sum = 0, sum2 = 0;
            for (double x : v) { sum += x; sum2 += x * x; }
            stats[m].mean = sum / n;
            stats[m].sd   = n > 1 ? std::sqrt((sum2 - sum*sum/n) / (n-1)) : 0.0;

            // Free raw data after finalisation
            v.clear();
            v.shrink_to_fit();
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Whole-structure metrics: cavities, patches, etc.
// Binned by protein size (total residues).
// ─────────────────────────────────────────────────────────────────────────────
enum WholeStructMetric
{
    WS_N_CAVITIES = 0,        // number of cavities
    WS_TOTAL_CAVITY_VOL,      // total cavity volume (ų)
    WS_LARGEST_CAVITY_VOL,    // largest single cavity (ų)
    WS_VOID_FRACTION,         // total_cavity / protein_volume

    WS_N_HPHOB_PATCHES,       // hydrophobic patch count (summed over chains)
    WS_LARGEST_HPHOB_SIZE,    // largest hydrophobic patch (residues)
    WS_LARGEST_HPHOB_KD,      // largest hydrophobic patch KD sum
    WS_N_POS_CLUSTERS,        // positive charge cluster count
    WS_LARGEST_POS_SIZE,      // largest positive cluster (residues)
    WS_LARGEST_POS_CHARGE,    // largest positive cluster net charge
    WS_N_NEG_CLUSTERS,        // negative charge cluster count
    WS_LARGEST_NEG_SIZE,      // largest negative cluster (residues)
    WS_LARGEST_NEG_CHARGE,    // largest negative cluster |net charge|
    WS_CLASHSCORE,            // heavy-atom clashes per 1000 atoms (MolProbity-style)

    WS_COUNT  // total number of whole-structure metrics
};

[[maybe_unused]] static const char* WholeStructMetricName(int m)
{
    static const char* names[] = {
        "n_cavities", "total_cavity_vol", "largest_cavity_vol", "void_fraction",
        "n_hphob_patches", "largest_hphob_size", "largest_hphob_kd",
        "n_pos_clusters", "largest_pos_size", "largest_pos_charge",
        "n_neg_clusters", "largest_neg_size", "largest_neg_charge",
        "clashscore"
    };
    return (m >= 0 && m < WS_COUNT) ? names[m] : "unknown";
}

// ── Protein size binning ────────────────────────────────────────────────────
static const int    N_SIZE_BINS = 5;
static const int    SIZE_EDGES[N_SIZE_BINS + 1] = {0, 100, 200, 400, 800, 999999};
//   bin 0: [0, 100)     small
//   bin 1: [100, 200)   medium-small
//   bin 2: [200, 400)   medium
//   bin 3: [400, 800)   large
//   bin 4: [800, ∞)     very large

static int SizeBin(int nResidues)
{
    for (int b = 0; b < N_SIZE_BINS; b++)
        if (nResidues < SIZE_EDGES[b + 1]) return b;
    return N_SIZE_BINS - 1;
}

static const char* SizeBinLabel(int b)
{
    static const char* labels[] = {
        "small(<100)", "med_small(100-200)", "medium(200-400)",
        "large(400-800)", "very_large(800+)"
    };
    return (b >= 0 && b < N_SIZE_BINS) ? labels[b] : "?";
}

// ── WholeStructCell: one size bin ──────────────────────────────────────────
struct WholeStructCell
{
    int count = 0;
    PercentileStats stats[WS_COUNT];
    std::vector<double> rawValues[WS_COUNT];

    void AddObservation(const double values[WS_COUNT])
    {
        count++;
        for (int m = 0; m < WS_COUNT; m++)
            if (!std::isnan(values[m]))
                rawValues[m].push_back(values[m]);
    }

    void Finalise()
    {
        static const double pctTargets[] = {0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};
        for (int m = 0; m < WS_COUNT; m++)
        {
            auto& v = rawValues[m];
            int n = (int)v.size();
            if (n == 0) continue;

            std::sort(v.begin(), v.end());
            for (int pi = 0; pi < 7; pi++)
            {
                double idx = pctTargets[pi] * (n - 1);
                int lo = (int)std::floor(idx);
                int hi = std::min(lo + 1, n - 1);
                double frac = idx - lo;
                stats[m].p[pi] = v[lo] * (1.0 - frac) + v[hi] * frac;
            }

            double sum = 0, sum2 = 0;
            for (double x : v) { sum += x; sum2 += x * x; }
            stats[m].mean = sum / n;
            stats[m].sd   = n > 1 ? std::sqrt((sum2 - sum*sum/n) / (n-1)) : 0.0;

            v.clear();
            v.shrink_to_fit();
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Ramachandran 2D grid, per amino acid type
//
// 72×72 bins covering phi ∈ [-180, +180) and psi ∈ [-180, +180) at 5° resolution.
// During the survey, raw counts are accumulated.  Finalise() converts counts
// to a density-based percentile for each bin:
//   percentile = (cumulative observations in cells of ≤ this density) / total × 100
//
// A residue in the most populated region scores ~100; an empty region scores 0.
// The quality checker then applies the standard surprise formula to this percentile.
// ─────────────────────────────────────────────────────────────────────────────
static const int    RAMA_N_BINS   = 72;        // per axis: 360° / 5° = 72
static const double RAMA_BIN_SIZE = 5.0;       // degrees

struct RamachandranGrid
{
    int      totalObs = 0;
    uint32_t counts[RAMA_N_BINS][RAMA_N_BINS] = {};      // accumulation phase
    double   pctile[RAMA_N_BINS][RAMA_N_BINS] = {};      // after Finalise()

    static int AngleBin(double angle)
    {
        // angle in [-180, +180); bin 0 covers [-180, -175)
        int b = (int)std::floor((angle + 180.0) / RAMA_BIN_SIZE);
        if (b < 0) b = 0;
        if (b >= RAMA_N_BINS) b = RAMA_N_BINS - 1;
        return b;
    }

    void AddObservation(double phi, double psi)
    {
        counts[AngleBin(phi)][AngleBin(psi)]++;
        totalObs++;
    }

    void Finalise()
    {
        if (totalObs == 0) return;

        // Collect (count, i, j) for all non-zero cells and sort by count ascending
        struct CellEntry { uint32_t cnt; int i, j; };
        std::vector<CellEntry> cells;
        cells.reserve(RAMA_N_BINS * RAMA_N_BINS);
        for (int i = 0; i < RAMA_N_BINS; i++)
            for (int j = 0; j < RAMA_N_BINS; j++)
                if (counts[i][j] > 0)
                    cells.push_back({counts[i][j], i, j});

        std::sort(cells.begin(), cells.end(),
                  [](const CellEntry& a, const CellEntry& b) { return a.cnt < b.cnt; });

        // Assign density-based percentile:
        // percentile = (cumulative observations in cells of ≤ this density) / total × 100
        int cumul = 0;
        for (auto& c : cells)
        {
            cumul += c.cnt;
            pctile[c.i][c.j] = 100.0 * (double)cumul / (double)totalObs;
        }
        // Cells with count=0 keep pctile=0 (completely empty = outlier)
    }

    // Return the density-based percentile for a (phi, psi) pair.
    // 0 = empty region (extreme outlier), ~100 = densest cluster.
    double Percentile(double phi, double psi) const
    {
        return pctile[AngleBin(phi)][AngleBin(psi)];
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// StructureReferenceDB: the full database
// ─────────────────────────────────────────────────────────────────────────────
class StructureReferenceDB
{
public:
    RefCell cells[N_AA][N_FACC_BINS];
    WholeStructCell wsCells[N_SIZE_BINS];   // whole-structure stats by protein size
    RamachandranGrid ramaGrids[N_AA];      // 2D (phi,psi) grids per AA
    int totalResidues = 0;
    int totalPDBs     = 0;
    std::mutex mtx_;   // for thread-safe AddObservation

    // ── Accumulation (survey phase) ──────────────────────────────────────

    // Thread-safe: adds one residue observation.
    // aa: one-letter code, facc: fractional accessibility [0,1]
    // values: array of SM_COUNT doubles (NaN for missing)
    void AddObservation(char aa, double facc, const double values[SM_COUNT])
    {
        int ai = AAIndex(aa);
        if (ai < 0) return;
        int bi = FaccBin(facc);

        std::lock_guard<std::mutex> lock(mtx_);
        cells[ai][bi].AddObservation(values);
        totalResidues++;
    }

    // Thread-safe: adds one (phi, psi) observation for the Ramachandran grid.
    void AddRamachandranObs(char aa, double phi, double psi)
    {
        int ai = AAIndex(aa);
        if (ai < 0) return;
        if (std::isnan(phi) || std::isnan(psi)) return;

        std::lock_guard<std::mutex> lock(mtx_);
        ramaGrids[ai].AddObservation(phi, psi);
    }

    // Thread-safe: adds one whole-structure observation.
    // nResidues: total residues in the structure (for size binning)
    // values: array of WS_COUNT doubles (NaN for missing)
    void AddWholeStructObs(int nResidues, const double values[WS_COUNT])
    {
        int si = SizeBin(nResidues);
        std::lock_guard<std::mutex> lock(mtx_);
        wsCells[si].AddObservation(values);
    }

    void Finalise()
    {
        for (int a = 0; a < N_AA; a++)
            for (int b = 0; b < N_FACC_BINS; b++)
                cells[a][b].Finalise();
        for (int s = 0; s < N_SIZE_BINS; s++)
            wsCells[s].Finalise();
        for (int a = 0; a < N_AA; a++)
            ramaGrids[a].Finalise();
    }

    // ── Lookup (checker phase) ───────────────────────────────────────────

    const RefCell* GetCell(char aa, double facc) const
    {
        int ai = AAIndex(aa);
        if (ai < 0) return nullptr;
        return &cells[ai][FaccBin(facc)];
    }

    const WholeStructCell* GetWSCell(int nResidues) const
    {
        return &wsCells[SizeBin(nResidues)];
    }

    // Get Ramachandran density-based percentile for a residue.
    // Returns -1 if AA is non-standard or grid has no data.
    double GetRamachandranPctile(char aa, double phi, double psi) const
    {
        int ai = AAIndex(aa);
        if (ai < 0) return -1.0;
        if (ramaGrids[ai].totalObs == 0) return -1.0;
        return ramaGrids[ai].Percentile(phi, psi);
    }

    bool HasRamachandranData() const
    {
        for (int a = 0; a < N_AA; a++)
            if (ramaGrids[a].totalObs > 0) return true;
        return false;
    }

    // ── Binary I/O ───────────────────────────────────────────────────────

    static const uint32_t MAGIC   = 0x42445253;  // "SRDB"
    static const uint32_t VERSION = 3;           // v3: added Ramachandran 2D grids

    bool WriteBinary(const std::string& path) const
    {
        std::ofstream f(path, std::ios::binary);
        if (!f) return false;

        uint32_t magic = MAGIC, ver = VERSION;
        uint32_t nAA = N_AA, nFacc = N_FACC_BINS, nMetrics = SM_COUNT;
        uint32_t nPdbs = (uint32_t)totalPDBs;
        uint32_t nRes  = (uint32_t)totalResidues;

        f.write((char*)&magic, 4);
        f.write((char*)&ver, 4);
        f.write((char*)&nAA, 4);
        f.write((char*)&nFacc, 4);
        f.write((char*)&nMetrics, 4);
        f.write((char*)&nPdbs, 4);
        f.write((char*)&nRes, 4);

        // Per-residue cells
        for (int a = 0; a < N_AA; a++)
        for (int b = 0; b < N_FACC_BINS; b++)
        {
            const RefCell& c = cells[a][b];
            uint32_t cnt = (uint32_t)c.count;
            f.write((char*)&cnt, 4);

            for (int m = 0; m < SM_COUNT; m++)
            {
                f.write((char*)c.stats[m].p, sizeof(double) * 7);
                f.write((char*)&c.stats[m].mean, sizeof(double));
                f.write((char*)&c.stats[m].sd, sizeof(double));
            }
        }

        // v2: whole-structure cells
        uint32_t nSizeBins = N_SIZE_BINS, nWSMetrics = WS_COUNT;
        f.write((char*)&nSizeBins, 4);
        f.write((char*)&nWSMetrics, 4);

        for (int s = 0; s < N_SIZE_BINS; s++)
        {
            const WholeStructCell& wc = wsCells[s];
            uint32_t cnt = (uint32_t)wc.count;
            f.write((char*)&cnt, 4);

            for (int m = 0; m < WS_COUNT; m++)
            {
                f.write((char*)wc.stats[m].p, sizeof(double) * 7);
                f.write((char*)&wc.stats[m].mean, sizeof(double));
                f.write((char*)&wc.stats[m].sd, sizeof(double));
            }
        }

        // v3: Ramachandran 2D grids
        uint32_t ramaN = RAMA_N_BINS;
        uint32_t ramaBinX100 = (uint32_t)(RAMA_BIN_SIZE * 100);  // 500
        f.write((char*)&ramaN, 4);
        f.write((char*)&ramaBinX100, 4);

        for (int a = 0; a < N_AA; a++)
        {
            const RamachandranGrid& rg = ramaGrids[a];
            uint32_t tot = (uint32_t)rg.totalObs;
            f.write((char*)&tot, 4);
            // Write the finalised percentile grid (72×72 doubles)
            f.write((char*)rg.pctile, sizeof(double) * RAMA_N_BINS * RAMA_N_BINS);
        }

        return f.good();
    }

    bool ReadBinary(const std::string& path)
    {
        std::ifstream f(path, std::ios::binary);
        if (!f) { std::cerr << "Cannot open " << path << std::endl; return false; }

        uint32_t magic, ver, nAA, nFacc, nMetrics, nPdbs, nRes;
        f.read((char*)&magic, 4);
        f.read((char*)&ver, 4);
        if (magic != MAGIC || ver < 1 || ver > VERSION)
        {
            std::cerr << "Invalid reference DB file (magic/version mismatch)" << std::endl;
            return false;
        }

        f.read((char*)&nAA, 4);
        f.read((char*)&nFacc, 4);
        f.read((char*)&nMetrics, 4);
        f.read((char*)&nPdbs, 4);
        f.read((char*)&nRes, 4);

        if (nAA != N_AA || nFacc != N_FACC_BINS || nMetrics != SM_COUNT)
        {
            std::cerr << "Reference DB dimension mismatch" << std::endl;
            return false;
        }

        totalPDBs     = (int)nPdbs;
        totalResidues = (int)nRes;

        // Per-residue cells
        for (int a = 0; a < N_AA; a++)
        for (int b = 0; b < N_FACC_BINS; b++)
        {
            RefCell& c = cells[a][b];
            uint32_t cnt;
            f.read((char*)&cnt, 4);
            c.count = (int)cnt;

            for (int m = 0; m < SM_COUNT; m++)
            {
                f.read((char*)c.stats[m].p, sizeof(double) * 7);
                f.read((char*)&c.stats[m].mean, sizeof(double));
                f.read((char*)&c.stats[m].sd, sizeof(double));
            }
        }

        // v2: whole-structure cells (absent in v1 files)
        if (ver >= 2)
        {
            uint32_t nSizeBins, nWSMetrics;
            f.read((char*)&nSizeBins, 4);
            f.read((char*)&nWSMetrics, 4);

            if (nSizeBins != N_SIZE_BINS || nWSMetrics != WS_COUNT)
            {
                std::cerr << "Reference DB whole-structure dimension mismatch" << std::endl;
                return false;
            }

            for (int s = 0; s < N_SIZE_BINS; s++)
            {
                WholeStructCell& wc = wsCells[s];
                uint32_t cnt;
                f.read((char*)&cnt, 4);
                wc.count = (int)cnt;

                for (int m = 0; m < WS_COUNT; m++)
                {
                    f.read((char*)wc.stats[m].p, sizeof(double) * 7);
                    f.read((char*)&wc.stats[m].mean, sizeof(double));
                    f.read((char*)&wc.stats[m].sd, sizeof(double));
                }
            }
        }
        else
        {
            std::cerr << "Note: v1 reference DB: no whole-structure distributions.\n";
        }

        // v3: Ramachandran 2D grids (absent in v1/v2 files)
        if (ver >= 3)
        {
            uint32_t ramaN, ramaBinX100;
            f.read((char*)&ramaN, 4);
            f.read((char*)&ramaBinX100, 4);

            if (ramaN != RAMA_N_BINS)
            {
                std::cerr << "Reference DB Ramachandran grid size mismatch ("
                          << ramaN << " vs " << RAMA_N_BINS << ")\n";
                return false;
            }

            for (int a = 0; a < N_AA; a++)
            {
                RamachandranGrid& rg = ramaGrids[a];
                uint32_t tot;
                f.read((char*)&tot, 4);
                rg.totalObs = (int)tot;
                f.read((char*)rg.pctile, sizeof(double) * RAMA_N_BINS * RAMA_N_BINS);
            }

            std::cerr << "Loaded Ramachandran grids (v3).\n";
        }
        else
        {
            std::cerr << "Note: v" << ver << " reference DB: no Ramachandran grids.\n";
        }

        return f.good();
    }

    // ── Human-readable summary ───────────────────────────────────────────

    void PrintSummary(std::ostream& out) const
    {
        out << "StructureReferenceDB: " << totalPDBs << " PDBs, "
            << totalResidues << " residues\n\n";

        out << "Observations per cell (AA × side-chain burial):\n";
        out << "AA";
        for (int b = 0; b < N_FACC_BINS; b++)
            out << "\t" << FaccBinLabel(b);
        out << "\ttotal\n";

        for (int a = 0; a < N_AA; a++)
        {
            out << AA_ORDER[a];
            int rowTotal = 0;
            for (int b = 0; b < N_FACC_BINS; b++)
            {
                out << "\t" << cells[a][b].count;
                rowTotal += cells[a][b].count;
            }
            out << "\t" << rowTotal << "\n";
        }

        out << "\nWhole-structure observations per size bin:\n";
        out << "Size_bin\tcount\n";
        for (int s = 0; s < N_SIZE_BINS; s++)
            out << SizeBinLabel(s) << "\t" << wsCells[s].count << "\n";

        out << "\nRamachandran observations per AA:\n";
        out << "AA\tobs\tnon_empty_cells\n";
        for (int a = 0; a < N_AA; a++)
        {
            int nonEmpty = 0;
            for (int i = 0; i < RAMA_N_BINS; i++)
                for (int j = 0; j < RAMA_N_BINS; j++)
                    if (ramaGrids[a].counts[i][j] > 0) nonEmpty++;
            out << AA_ORDER[a] << "\t" << ramaGrids[a].totalObs
                << "\t" << nonEmpty << "/" << (RAMA_N_BINS * RAMA_N_BINS) << "\n";
        }
    }
};
