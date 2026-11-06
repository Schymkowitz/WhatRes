// SequenceAnalysis, Copyright Switch Lab 2017
// Joost Schymkowitz

#ifndef PDB_FILE_H
#define PDB_FILE_H

#include <sstream>
#include <iostream>
#include <fstream>

#include "Peptide.h"
using namespace std;


class PDB_File
{
public:
	vector<Peptide> peptides;
	vector<string> sequences;
	vector<char> chain_names;
	string pdbName;
	int current_P;

	// Header records preserved verbatim for downstream tools that need
	// to parse them (e.g. BiomtExpander reads remark350Lines for BIOMT
	// operators; CrystalLayerExpander reads remark290Lines and cryst1Line).
	// PDB_Reader populates these during loading.
	vector<string> remark350Lines;   // BIOMT (biological assembly)
	vector<string> remark290Lines;   // SMTRY (crystallographic symmetry)
	string         cryst1Line;       // single CRYST1 record
	int            nModels;          // count of MODEL records (NMR signal)

	// Experimental metadata captured from header records.
	string         experimentalMethod;  // EXPDTA, e.g. "X-RAY DIFFRACTION"
	double         resolutionA;         // REMARK 2 resolution (Å); -1 if unknown

	// Provenance lines (REMARK 999) accumulated by pipeline stages
	// (BiomtExpander, CrystalLayerExpander, LayerNormalizer, ingest_structure_command).
	// Each entry is the body text after "REMARK 999 ", emitted by
	// WritePDBFormat as full REMARK 999 records at the top of the output.
	vector<string> provenanceLines;

	PDB_File(): current_P(0), nModels(0), resolutionA(-1.0)
	  {}
  	
	void Add_Peptide(Peptide& P) 	{ peptides.push_back(P); };
	void Set_currentP(int x)	{ current_P = x; }		
	bool Set_currentP(char a);
	void Set_lastP_as_current()	{ current_P = peptides.size()-1; }		
	int Get_currentP() 	{ return current_P; }
	string Get_currentSeq() { return sequences[current_P];}
	void Add_Res_to_currentP(Residue& R) { peptides[current_P].Add_Residue(R); }
	int Pepcount() { return peptides.size()-1; }
	int CurrentP_ResCount() { return peptides[current_P].Rescount(); }
	void Add_Atom_to_lastrest(Atom& A) { peptides[current_P].Add_Atom_to_lastrest(A); }
	void Sequence_on_screen();
	int identify_chain(const char c);
	int Get_Pos_from_PDBnumber(int pdbnum, const char c);
	string GetChainSeq(const char c);
	string GetChainSeq(int pepnum);
	double dist(Atom& a1, Atom& a2);
	void SetContactOrder(double dist_cutoff, double seqdist_cutoff);
	void PrintDihedrals(ofstream& out);
	void WritePDBFormat(ofstream& out);
	void WriteChainPDBFormat(ofstream& out, char c);
	bool CheckChain(char c);
};

inline ostream& operator<<(ostream& os, const PDB_File& pdb)  
{  
    os << pdb.pdbName << endl;
	vector<Peptide>::const_iterator it = pdb.peptides.begin(), ite = pdb.peptides.end();
	for(; it!=ite; it++) os << (*it);
    return os;  
} 

#endif
