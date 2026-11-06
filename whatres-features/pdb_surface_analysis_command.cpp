// pdb_surface_analysis_command.cpp
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz

#include "pdb_surface_analysis_command.h"
#include "SurfacePatcher.h"
#include "CavityDetector.h"
#include "InterfaceVoidAnalyser.h"
#include "PDB_File.h"
#include "PDB_Reader.h"
#include "LineReader.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

// ═════════════════════════════════════════════════════════════════════════════
// Surface patch output helpers (unchanged)
// ═════════════════════════════════════════════════════════════════════════════

static void writeHeader(ofstream& out)
{
    out << "pdb_name\tchain"         // file + chain identity
        << "\tn_res"                  // total ATOM residues in chain
        << "\tn_surf"                 // exposed residues
        << "\tfrac_surf"             // fraction exposed
        << "\tmean_burial"           // mean Cβ-neighbour count (all residues)
        << "\tsurf_net_charge"       // net charge of exposed residues
        << "\tfrac_surf_hphob"       // fraction of exposed residues that are hydrophobic
        // ── Hydrophobic patches ──────────────────────────────────────────
        << "\tn_hphob_patches"       // number of patches >= min_patch_size
        << "\thphob1_size"           // residues in largest patch
        << "\thphob1_kd_sum"         // Kyte-Doolittle sum of largest patch
        << "\thphob1_rg"             // radius of gyration (Å)
        << "\thphob1_aromatic"       // aromatic (F,W,Y,H) count
        << "\thphob1_entropy"        // Shannon composition entropy (bits)
        << "\thphob1_mean_burial"    // mean burial depth of patch members
        // ── Positive charge clusters ─────────────────────────────────────
        << "\tn_pos_clusters"        // number of positive clusters
        << "\tpos1_size"
        << "\tpos1_charge"           // net charge of largest positive cluster
        << "\tpos1_rg"
        << "\tpos1_mean_burial"
        // ── Negative charge clusters ─────────────────────────────────────
        << "\tn_neg_clusters"
        << "\tneg1_size"
        << "\tneg1_charge"
        << "\tneg1_rg"
        << "\tneg1_mean_burial"
        // ── Residue lists ────────────────────────────────────────────────
        << "\thphob1_residues"       // e.g. "I42,L45,V48"
        << "\tpos1_residues"
        << "\tneg1_residues"
        << "\n";
}

static void writePatchFields(ofstream& out, const SurfacePatch& p, int patchCount)
{
    out << "\t" << patchCount
        << "\t" << p.size
        << "\t" << fixed << setprecision(2) << p.kd_sum
        << "\t" << fixed << setprecision(2) << p.rg
        << "\t" << p.aromatic
        << "\t" << fixed << setprecision(3) << p.comp_entropy
        << "\t" << fixed << setprecision(1) << p.mean_burial;
}

static void writeChargeClusterFields(ofstream& out, const SurfacePatch& p, int clusterCount)
{
    out << "\t" << clusterCount
        << "\t" << p.size
        << "\t" << fixed << setprecision(2) << p.net_charge
        << "\t" << fixed << setprecision(2) << p.rg
        << "\t" << fixed << setprecision(1) << p.mean_burial;
}

static void writeRow(ofstream& out, const SurfaceAnalysis& r)
{
    out << r.pdb_name
        << "\t" << r.chain
        << "\t" << r.n_total
        << "\t" << r.n_surface
        << "\t" << fixed << setprecision(3) << r.frac_surface
        << "\t" << fixed << setprecision(1) << r.mean_burial
        << "\t" << fixed << setprecision(2) << r.surf_net_charge
        << "\t" << fixed << setprecision(3) << r.frac_surf_hphob;

    writePatchFields(out, r.topHydrophobic(), (int)r.hydrophobic.size());
    writeChargeClusterFields(out, r.topPositive(), (int)r.positive.size());
    writeChargeClusterFields(out, r.topNegative(), (int)r.negative.size());

    out << "\t" << r.topHydrophobic().residue_list
        << "\t" << r.topPositive().residue_list
        << "\t" << r.topNegative().residue_list
        << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Verbose per-patch detail file
// ─────────────────────────────────────────────────────────────────────────────
static void writePatchDetailHeader(ofstream& out)
{
    out << "pdb_name\tchain\tpatch_type\trank\tsize"
        << "\tkd_sum\tnet_charge\trg\taromatic\tentropy\tmean_burial\tresidues\n";
}

static void writePatchDetailRows(ofstream& out, const SurfaceAnalysis& r,
                                  const string& type, const vector<SurfacePatch>& patches)
{
    for (int i = 0; i < (int)patches.size(); i++)
    {
        const SurfacePatch& p = patches[i];
        out << r.pdb_name << "\t" << r.chain << "\t" << type << "\t" << (i + 1)
            << "\t" << p.size
            << "\t" << fixed << setprecision(2) << p.kd_sum
            << "\t" << fixed << setprecision(2) << p.net_charge
            << "\t" << fixed << setprecision(2) << p.rg
            << "\t" << p.aromatic
            << "\t" << fixed << setprecision(3) << p.comp_entropy
            << "\t" << fixed << setprecision(1) << p.mean_burial
            << "\t" << p.residue_list
            << "\n";
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Cavity output helpers
// ═════════════════════════════════════════════════════════════════════════════

static void writeCavitySummaryHeader(ofstream& out)
{
    out << "pdb_name"
        << "\tn_cavities"
        << "\ttotal_cavity_vol"
        << "\tlargest_cavity_vol"
        << "\tprotein_volume"
        << "\tvoid_fraction"
        << "\n";
}

static void writeCavitySummaryRow(ofstream& out, const string& pdbName,
                                   const CavitySummary& cs)
{
    out << pdbName
        << "\t" << cs.n_cavities
        << "\t" << fixed << setprecision(1) << cs.total_volume
        << "\t" << fixed << setprecision(1) << cs.largest_volume
        << "\t" << fixed << setprecision(1) << cs.protein_volume
        << "\t" << fixed << setprecision(4) << cs.void_fraction
        << "\n";
}

static void writeCavityDetailHeader(ofstream& out)
{
    out << "pdb_name\tcavity_id\tvolume\tcentroid_x\tcentroid_y\tcentroid_z"
        << "\tn_lining_residues\tlining_residues\n";
}

static void writeCavityDetailRows(ofstream& out, const string& pdbName,
                                    const CavitySummary& cs)
{
    for (const auto& c : cs.cavities)
    {
        out << pdbName
            << "\t" << c.id
            << "\t" << fixed << setprecision(1) << c.volume
            << "\t" << fixed << setprecision(2) << c.cx
            << "\t" << fixed << setprecision(2) << c.cy
            << "\t" << fixed << setprecision(2) << c.cz
            << "\t" << (int)c.lining_residues.size();

        // Residue list
        string resList;
        for (const auto& r : c.lining_residues)
        {
            if (!resList.empty()) resList += ',';
            resList += r;
        }
        out << "\t" << resList << "\n";
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Interface void output helpers
// ═════════════════════════════════════════════════════════════════════════════

static void writeInterfaceHeader(ofstream& out)
{
    out << "pdb_name"
        << "\tchainA\tchainB"
        << "\tn_iface_res_A\tn_iface_res_B"
        << "\tiface_void_vol"
        << "\tiface_total_vol"
        << "\tiface_occupied_vol"
        << "\tpacking_density"
        << "\tvoid_index"
        << "\tiface_residues_A"
        << "\tiface_residues_B"
        << "\n";
}

static void writeInterfaceRow(ofstream& out, const string& pdbName,
                                const InterfaceVoidResult& r)
{
    out << pdbName
        << "\t" << r.chainA
        << "\t" << r.chainB
        << "\t" << r.n_interface_A
        << "\t" << r.n_interface_B
        << "\t" << fixed << setprecision(1) << r.interface_void_volume
        << "\t" << fixed << setprecision(1) << r.interface_total_volume
        << "\t" << fixed << setprecision(1) << r.interface_occupied
        << "\t" << fixed << setprecision(4) << r.packing_density
        << "\t" << fixed << setprecision(2) << r.void_index
        << "\t" << r.residues_A
        << "\t" << r.residues_B
        << "\n";
}

// ═════════════════════════════════════════════════════════════════════════════
// Execute
// ═════════════════════════════════════════════════════════════════════════════
void pdb_surface_analysis_command::Execute()
{
    if (Option.inputfiles.empty())
    {
        cout << "pdb surface analysis: no input files specified." << endl;
        return;
    }

    // ── Configure engines ────────────────────────────────────────────────
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
    cavDet.lining_dist  = Option.cavity_lining_dist;

    InterfaceVoidAnalyser ifaceAnalyser;
    ifaceAnalyser.grid_spacing  = Option.interface_grid_spacing;
    ifaceAnalyser.contact_dist  = Option.interface_contact_dist;
    ifaceAnalyser.shell_dist    = Option.interface_shell_dist;
    ifaceAnalyser.probe_radius  = Option.cavity_probe_radius;   // reuse probe
    ifaceAnalyser.box_padding   = Option.interface_box_padding;

    // ── Open output files ────────────────────────────────────────────────
    // 1. Surface patches (per-chain)
    string surfPath = Option.prefix + "_surface_analysis.txt";
    ofstream surfOut(surfPath.c_str());
    writeHeader(surfOut);

    // 2. Per-patch detail (optional)
    ofstream detailOut;
    if (Option.surface_verbose)
    {
        string detailPath = Option.prefix + "_surface_patches.txt";
        detailOut.open(detailPath.c_str());
        writePatchDetailHeader(detailOut);
    }

    // 3. Cavity summary (per-PDB)
    string cavSumPath = Option.prefix + "_cavities.txt";
    ofstream cavSumOut(cavSumPath.c_str());
    writeCavitySummaryHeader(cavSumOut);

    // 4. Cavity detail (per-cavity, optional, written when verbose)
    ofstream cavDetOut;
    if (Option.surface_verbose)
    {
        string cavDetPath = Option.prefix + "_cavity_details.txt";
        cavDetOut.open(cavDetPath.c_str());
        writeCavityDetailHeader(cavDetOut);
    }

    // 5. Interface void (per-chain-pair)
    string ifacePath = Option.prefix + "_interface_voids.txt";
    ofstream ifaceOut(ifacePath.c_str());
    writeInterfaceHeader(ifaceOut);

    // ── Loop over PDB list ───────────────────────────────────────────────
    int nPdbs = 0;
    for (const string& listFile : Option.inputfiles)
    {
        ifstream list(listFile.c_str());
        if (!list)
        {
            cout << "pdb surface analysis: cannot open list file: " << listFile << endl;
            continue;
        }

        string pdbPath;
        while (getline(list, pdbPath))
        {
            if (!pdbPath.empty() && pdbPath.back() == '\r')
                pdbPath.pop_back();
            if (pdbPath.empty()) continue;

            // Prepend pdb directory for relative paths
            if (!Option.surface_pdb_dir.empty() && pdbPath[0] != '/' && pdbPath[0] != '~')
            {
                string dir = Option.surface_pdb_dir;
                if (dir.back() != '/') dir += '/';
                pdbPath = dir + pdbPath;
            }

            PDB_File pdb;
            PDB_Reader reader;
            if (reader.Readpdb(pdbPath, pdb) < 0)
            {
                cout << "pdb surface analysis: could not read " << pdbPath << ", skipped." << endl;
                continue;
            }

            // Strip directory prefix for display name
            string pdbName = pdbPath;
            size_t slash = pdbName.find_last_of("/\\");
            if (slash != string::npos) pdbName = pdbName.substr(slash + 1);

            nPdbs++;
            cerr << " " << pdbName << flush;

            // ── 1. Surface patches (per-chain) ──────────────────────────
            for (int c = 0; c < (int)pdb.peptides.size(); c++)
            {
                if (pdb.peptides[c].reslist.empty()) continue;
                SurfaceAnalysis result = patcher.Analyse(pdb.peptides[c], pdbName);
                writeRow(surfOut, result);

                if (Option.surface_verbose)
                {
                    writePatchDetailRows(detailOut, result, "HYDROPHOBIC", result.hydrophobic);
                    writePatchDetailRows(detailOut, result, "POSITIVE",    result.positive);
                    writePatchDetailRows(detailOut, result, "NEGATIVE",    result.negative);
                }
            }

            // ── 2. Cavity detection (whole PDB) ─────────────────────────
            CavitySummary cavResult = cavDet.Detect(pdb);
            writeCavitySummaryRow(cavSumOut, pdbName, cavResult);

            if (Option.surface_verbose)
                writeCavityDetailRows(cavDetOut, pdbName, cavResult);

            // ── 3. Interface void analysis (all chain pairs) ────────────
            vector<InterfaceVoidResult> ifaceResults = ifaceAnalyser.AnalyseAll(pdb);
            for (const auto& ir : ifaceResults)
                writeInterfaceRow(ifaceOut, pdbName, ir);
        }
    }

    cerr << endl;
    cout << "pdb surface analysis: " << nPdbs << " PDBs processed." << endl;
    cout << "  surface patches:  " << surfPath << endl;
    cout << "  cavities:         " << cavSumPath << endl;
    cout << "  interface voids:  " << ifacePath << endl;
    if (Option.surface_verbose)
        cout << "  (verbose details also written)" << endl;
}
