// optionclass.cpp
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz

#include "optionclass.h"
#include <cstdlib>
#include <fstream>
#include <vector>


void Option_class::LoadDefaultOptionSet()
{
	pH = 7.4;
	inputfiles.clear();
	mutations.clear();
	prefix = "PSX";

	// PDB surface analysis
	surface_burial_radius    = 10.0;
	surface_burial_threshold = 16;
	surface_adj_cutoff       = 8.0;
	surface_min_patch_size   = 2;
	surface_hphob_set        = "AVILMFWYPCG";
	surface_verbose          = false;
	surface_pdb_dir          = "";

	// Cavity detection
	cavity_grid_spacing      = 1.0;
	cavity_probe_radius      = 1.4;
	cavity_min_volume        = 10.0;
	cavity_lining_dist       = 3.5;

	// Interface void analysis
	interface_grid_spacing   = 0.5;
	interface_contact_dist   = 5.0;
	interface_shell_dist     = 4.0;
	interface_box_padding    = 3.0;

	// PDB quality check
	qc_reference_db      = "";
	qc_outlier_lo_pctile = 0;   // p1
	qc_outlier_hi_pctile = 6;   // p99
	qc_min_ref_count     = 50;
}

bool Option_class::Read_Options(string filename)
{
	OptionFile X;
	if (!X.Read(filename))
	{
		cout<<"Optionfile not found, using defaults"<<endl;
		return false;
	}
	string line;
	X.Get("TITLE", line);
	// WHATRES_optionfile; AGADIR_optionfile is accepted for older options files
	if (line != "WHATRES_optionfile" && line != "AGADIR_optionfile")
	{
		cout << "Not a valid optionfile - defaults kept." << endl;
		return false;
	}

	X.Get("pH", pH);
	X.GetAll("inputfiles", inputfiles);
	X.GetAll("mutations", mutations);
	X.Get("prefix", prefix);

	// PDB surface analysis
	X.Get("surface_burial_radius", surface_burial_radius);
	X.Get("surface_burial_threshold", surface_burial_threshold);
	X.Get("surface_adj_cutoff", surface_adj_cutoff);
	X.Get("surface_min_patch_size", surface_min_patch_size);
	X.Get("surface_hphob_set", surface_hphob_set);
	X.Get("surface_verbose", surface_verbose);
	X.Get("surface_pdb_dir", surface_pdb_dir);

	// Cavity detection
	X.Get("cavity_grid_spacing", cavity_grid_spacing);
	X.Get("cavity_probe_radius", cavity_probe_radius);
	X.Get("cavity_min_volume", cavity_min_volume);
	X.Get("cavity_lining_dist", cavity_lining_dist);

	// Interface void analysis
	X.Get("interface_grid_spacing", interface_grid_spacing);
	X.Get("interface_contact_dist", interface_contact_dist);
	X.Get("interface_shell_dist", interface_shell_dist);
	X.Get("interface_box_padding", interface_box_padding);

	// PDB quality check
	X.Get("qc_reference_db", qc_reference_db);
	X.Get("qc_outlier_lo_pctile", qc_outlier_lo_pctile);
	X.Get("qc_outlier_hi_pctile", qc_outlier_hi_pctile);
	X.Get("qc_min_ref_count", qc_min_ref_count);

	return true;
}

void Option_class::ReplaceInputFiles(vector<string>& newfiles)
{
	inputfiles.clear();
	for(int i=0;i<newfiles.size();i++) inputfiles.push_back(newfiles[i]);
}

void Option_class::ReplaceInputFiles(const string& newfile)
{
	inputfiles.clear();
	inputfiles.push_back(newfile);
}

void Option_class::ReplaceMutations(vector<string>& newmuts)
{
	mutations.clear();
	for(int i=0;i<newmuts.size();i++) mutations.push_back(newmuts[i]);
}

void Option_class::ReplacePrefix(const string& pre)
{
	prefix = pre;
}

void Option_class::operator=(const Option_class &O)
{
	inputfiles = O.inputfiles;
	mutations  = O.mutations;
	prefix     = O.prefix;
	pH         = O.pH;

	surface_burial_radius    = O.surface_burial_radius;
	surface_burial_threshold = O.surface_burial_threshold;
	surface_adj_cutoff       = O.surface_adj_cutoff;
	surface_min_patch_size   = O.surface_min_patch_size;
	surface_hphob_set        = O.surface_hphob_set;
	surface_verbose          = O.surface_verbose;
	surface_pdb_dir          = O.surface_pdb_dir;

	cavity_grid_spacing      = O.cavity_grid_spacing;
	cavity_probe_radius      = O.cavity_probe_radius;
	cavity_min_volume        = O.cavity_min_volume;
	cavity_lining_dist       = O.cavity_lining_dist;

	interface_grid_spacing   = O.interface_grid_spacing;
	interface_contact_dist   = O.interface_contact_dist;
	interface_shell_dist     = O.interface_shell_dist;
	interface_box_padding    = O.interface_box_padding;

	qc_reference_db      = O.qc_reference_db;
	qc_outlier_lo_pctile = O.qc_outlier_lo_pctile;
	qc_outlier_hi_pctile = O.qc_outlier_hi_pctile;
	qc_min_ref_count     = O.qc_min_ref_count;
}

Option_class::Option_class(const Option_class &O)
{
	*this = O;
}

void Option_class::PrintOptions(ofstream &o)
{
	o << PrintOptions().str();
}

stringstream Option_class::PrintOptions()
{
	stringstream o;
	o << "pH:" << "\t" << pH << endl;
	o << "inputfiles: " << "\t";
	for(int i=0;i<inputfiles.size();i++) o << inputfiles[i] << "\t";
	o << "\n";
	o << "prefix:" << "\t" << prefix << endl;
	o << "surface_burial_radius:"    << "\t" << surface_burial_radius    << endl;
	o << "surface_burial_threshold:" << "\t" << surface_burial_threshold << endl;
	o << "surface_adj_cutoff:"       << "\t" << surface_adj_cutoff       << endl;
	o << "surface_min_patch_size:"   << "\t" << surface_min_patch_size   << endl;
	o << "surface_hphob_set:"        << "\t" << surface_hphob_set        << endl;
	o << "surface_verbose:"          << "\t" << surface_verbose          << endl;
	o << "surface_pdb_dir:"          << "\t" << surface_pdb_dir          << endl;
	o << "cavity_grid_spacing:"      << "\t" << cavity_grid_spacing      << endl;
	o << "cavity_probe_radius:"      << "\t" << cavity_probe_radius      << endl;
	o << "cavity_min_volume:"        << "\t" << cavity_min_volume        << endl;
	o << "cavity_lining_dist:"       << "\t" << cavity_lining_dist       << endl;
	o << "interface_grid_spacing:"   << "\t" << interface_grid_spacing   << endl;
	o << "interface_contact_dist:"   << "\t" << interface_contact_dist   << endl;
	o << "interface_shell_dist:"     << "\t" << interface_shell_dist     << endl;
	o << "interface_box_padding:"    << "\t" << interface_box_padding    << endl;
	o << "qc_reference_db:"          << "\t" << qc_reference_db          << endl;
	o << "qc_outlier_lo_pctile:"     << "\t" << qc_outlier_lo_pctile     << endl;
	o << "qc_outlier_hi_pctile:"     << "\t" << qc_outlier_hi_pctile     << endl;
	o << "qc_min_ref_count:"         << "\t" << qc_min_ref_count         << endl;
	return o;
}
