// SequenceAnalysis, Copyright Switch Lab 2012
// Joost Schymkowitz

#ifndef PDB_READER_H
#define PDB_READER_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include "PDB_LineInterpretation.h"
#include "PDB_File.h"

using namespace std;

class PDB_Reader
{
public:
	ifstream pdbfile;
	PDB_LineInterpretation LIP;
	int MolNum;
	stringstream unknown;
	string pdbName;

	PDB_Reader() { MolNum = 0; }

	int OpenPdb(const string &naam)
	{
		pdbName=naam;
		pdbfile.open(naam.c_str(), ios::binary);
		if (!pdbfile) { 
			cerr<<"can't open pdb-file: " <<naam<<endl; 
			return -1;
		}
		return 1;
	}

	void ClosePdb() { pdbfile.close(); }

	int Readpdb(const string &naam, PDB_File& pdb);
	
};


#endif
