// optionclass.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Configuration object for the structure-quality commands: the options read
// by pdb_quality_check, pdb_surface_analysis, PDB_dist and structure_survey.
//
// The options file is a list of <tag>value; entries, starting with
// <TITLE>WHATRES_optionfile; and ending with <END>. Fields absent from the file
// keep their defaults (see LoadDefaultOptionSet). Tags this class does not
// know are ignored.
#ifndef OPTION_CLASS_H
#define OPTION_CLASS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "OptionFile.h"

using namespace std;

class Option_class
{
public:
	// I/O
	vector<string> inputfiles;  // input files (from -in or <inputfiles> tags); here a PDB list file
	vector<string> mutations;   // set by -mut / <mutations>; not used by these commands
	string prefix;              // output file name prefix (default "PSX"; override with -pre)

	// Conditions
	double pH;                  // solution pH (default 7.4); charge states in SurfacePatcher

	// PDB surface analysis: SurfacePatcher
	double surface_burial_radius;    // CB-neighbour sphere radius (A); residues with fewer
	                                 // than surface_burial_threshold neighbours are "exposed"
	int    surface_burial_threshold; // exposure cutoff: exposed if neighbour count < this
	double surface_adj_cutoff;       // CB-CB distance (A) for patch adjacency graph edges
	int    surface_min_patch_size;   // minimum connected residues to report a patch
	string surface_hphob_set;        // one-letter codes treated as hydrophobic
	bool   surface_verbose;          // if true: also write per-patch detail file (*_surface_patches.txt)
	string surface_pdb_dir;          // directory prepended to relative paths in the PDB list file (optional)

	// Cavity detection
	double cavity_grid_spacing;      // voxel edge length for cavity detection (A); smaller = finer
	double cavity_probe_radius;      // solvent probe radius for volume exclusion (A)
	double cavity_min_volume;        // minimum cavity volume to report (A^3)
	double cavity_lining_dist;       // max atom-centroid distance for lining residues (A)

	// Interface void analysis
	double interface_grid_spacing;   // voxel edge for interface box (A)
	double interface_contact_dist;   // cutoff to define interface residues (A)
	double interface_shell_dist;     // void voxel must have both chains within this (A)
	double interface_box_padding;    // padding around interface residues (A)

	// PDB quality check
	string qc_reference_db;          // path to structure_reference.bin
	int    qc_outlier_lo_pctile;     // lower outlier percentile index (0=p1, 1=p5, 2=p25)
	int    qc_outlier_hi_pctile;     // upper outlier percentile index (4=p75, 5=p95, 6=p99)
	int    qc_min_ref_count;         // minimum observations in ref cell to trust

	Option_class() { LoadDefaultOptionSet(); }
	Option_class(const Option_class &O);   // copies every field
	~Option_class() {}

	void LoadDefaultOptionSet();
	bool Read_Options(string filename);    // false if the file is missing or not an options file (defaults kept)
	void operator=(const Option_class &O);

	void PrintInputFiles()
	{
		cout << "inputfiles: ";
		for(int i = 0; i < (int)inputfiles.size(); i++) cout << inputfiles[i] << " ";
		cout << endl;
	}
	void ReplaceInputFiles(vector<string>& newfiles);
	void ReplaceInputFiles(const string& newfile);
	void ReplaceMutations(vector<string>& newmuts);
	void ReplacePrefix(const string& pre);

	void PrintOptions(ofstream &o);        // writes PSX_usedOptions.txt
	stringstream PrintOptions();
};

#endif
