// SequenceAnalysis, Copyright Switch Lab 2012
// Joost Schymkowitz

#include "PDB_LineInterpretation.h"

bool PDB_LineInterpretation::standardise_atom()
{
	
	// Look up the names of atoms to verify if they conform to our standard
	string key = amino + atom;
	string temp_amino;
//	tablesisolvation *test;
	if (amino == "CSE")
	{
		amino = "CYS"; // seleno cysteine ->convert to ATOM, CYS
		keyword = "ATOM";
		if (atom=="E" || atom=="SE")
			atom = "SG";
	}
	else if (amino == "MSE")
	{
		amino = "MET"; // seleno methionine -> convert to ATOM, MET
		keyword = "ATOM";
		if (atom=="E" || atom=="SE") 
			atom = "SD";
	}
	else if (amino == "GNP")
	{
		amino = "GTP"; // GMPNP -> convert to GTP
		if (atom=="N3B") 
			atom = "O3B";
	}
	else if (amino=="DG")
	{
		amino="G";
	}
	else if (amino=="DC")
	{
		amino="C";
	}
	else if (amino=="DT")
	{
		amino="T";
	}
	else if (amino=="DA")
	{
		amino="A";
	}
	if (atom=="C7" && amino=="T")
	{atom="C5M";    }
	temp_amino = amino;
	if (atom == "OT1")
		atom = "O";
	if (atom == "OT2")
		atom = "OXT";
	if (amino =="HOH" && keyword == "ATOM")
		keyword = "HETATM"; 
	if (amino == "H2O")
		amino = "HOH";
	/*
    atom = LOOKUP.GetStandardAtom("XXX", atom);
	atom = LOOKUP.GetStandardAtom(amino, atom);
	residue_type resType=LOOKUP.GetResidueType(amino);
	if (resType != SINGLEATOMLIGAND && resType != UNDEFINED && resType != WATER)
		keyword="ATOM";

	if ((resType==PROTEIN && (atom == "OXT" || atom == "tN")) || (resType==DNA && atom == "O3P")) 
		temp_amino = "XXX";
	//cout << amino<<" "<<atom<<endl;
	test = LOOKUP.GetSolvation(temp_amino, atom);
	if (test==NULL && !proton)
	{
		ofstream x ("Unrecognized_molecules.txt", ios::app);
		//cerr << "Unrecognised atom: " << temp_amino <<" "<< resnumber << " " << atom << endl;
		x<< "Unrecognised atom: " << temp_amino <<" "<< resnumber << " " << atom << endl;
		x.close();
		return false;
	}
	if (proton)
		return false; */
	return true;
}
