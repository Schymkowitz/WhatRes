// SequenceAnalysis, Copyright Switch Lab 2017
// Joost Schymkowitz

#ifndef RESIDUE_H
#define RESIDUE_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include "Atom.h"
#include "vector3D.h"

using namespace std;

class Residue
{
public:
	string molecule_type;
	string residue_name;
	char chain;
	vector<Atom> atoms;
	int C, N, CA; 	// index to backbone atoms for dehidrals
	int number_in_pdb;
	double ContactOrder;
	double phi;
	double psi;
	double omega;
	vector<int> seqdist; // list of absolute sequence distance of its interactors
	vector<int> abs_int; // list of sequence position of the interactors
	int max_interactor;	// the interactor that has the highest position in the chain
	int interactor_mean; // mean sequence distance interactor
	int interactor_median; // median sequence distance interactor
    int interactor_80; // interactor that completes 80% of iteractions
	
	
	Residue(const string& MT, const string& RN, char C, int NIP)
		: molecule_type(MT), residue_name(RN), chain(C),
		  C(-1), N(-1), CA(-1),
		  number_in_pdb(NIP), ContactOrder(0.),
		  phi(-999.), psi(-999.), omega(-999.),
		  max_interactor(-1)
	{
	}
		
	void Add_Atom(Atom& A)
	{
		int i = atoms.size();
		atoms.push_back(A);
		if(A.atom_name == "N") N = i;
		else if(A.atom_name == "CA") CA = i;
		else if(A.atom_name == "C") C = i;
	}
	
	string ResID()
	{
		stringstream s;
		s << residue_name << chain << number_in_pdb;
		return s.str();
	}
	
	void TruncateToAlanine()
	{
	    // Already minimal
	    if(residue_name == "ALA" || residue_name == "GLY")
	        return;

	    std::vector<Atom> kept;
	    kept.reserve(5);

	    for(auto &a : atoms)
	    {
	        if(a.atom_name == "N"  ||
	           a.atom_name == "CA" ||
	           a.atom_name == "C"  ||
	           a.atom_name == "O"  ||
	           a.atom_name == "CB")
	        {
	            Atom aa = a;
	            kept.push_back(aa);
	        }
	    }

	    atoms.swap(kept);

	    // update residue name at residue level
	    residue_name = "ALA";
	}
	
	char OneLetter() const
	{
	    // residue_name is expected to be 3-letter (e.g. "ALA")
	    const std::string& r = residue_name;

	    if(r == "ALA") return 'A';
	    if(r == "ARG") return 'R';
	    if(r == "ASN") return 'N';
	    if(r == "ASP") return 'D';
	    if(r == "CYS") return 'C';
	    if(r == "GLN") return 'Q';
	    if(r == "GLU") return 'E';
	    if(r == "GLY") return 'G';
	    if(r == "HIS") return 'H';
	    if(r == "ILE") return 'I';
	    if(r == "LEU") return 'L';
	    if(r == "LYS") return 'K';
	    if(r == "MET") return 'M';
	    if(r == "PHE") return 'F';
	    if(r == "PRO") return 'P';
	    if(r == "SER") return 'S';
	    if(r == "THR") return 'T';
	    if(r == "TRP") return 'W';
	    if(r == "TYR") return 'Y';
	    if(r == "VAL") return 'V';

	    // common alternates you might encounter
	    if(r == "MSE") return 'M'; // selenomethionine
	    if(r == "SEC") return 'U'; // selenocysteine
	    if(r == "PYL") return 'O'; // pyrrolysine

	    return 'X';
	}
	Atom& getC() { return atoms[C];}
	Atom& getCA() { return atoms[CA];}
	Atom& getN() { return atoms[N];}
	bool checkC() { return (C != -1);}
	bool checkCA() { return (CA != -1);}
	bool checkN() { return (N != -1);}
	
	void WritetoPDB(ofstream& o)
	{
		int sizer = atoms.size();
		for(int i=0; i < sizer; i++)
		{
			atoms[i].PutPDBline(o, 0., number_in_pdb, i, chain, residue_name);
		}
	}
	
	void WritetoPDB(std::ostream& o)
	{
	    int sizer = atoms.size();
	    for(int i=0; i < sizer; i++)
	    {
	        atoms[i].PutPDBline(o, 0., number_in_pdb, i, chain, residue_name);
	    }
	}
};

inline ostream& operator<<(ostream& os, const Residue& res)  
{  
    os << res.molecule_type << "\t" << res.residue_name << "\t" << res.chain << "\t" << res.number_in_pdb << endl;
    vector<Atom>::const_iterator it = res.atoms.begin(), ite = res.atoms.end();
	for(; it!=ite; it++) os << (*it);
    return os;  
} 

inline double Dihedral(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4)
{
	Vector3D i(a1.x, a1.y, a1.z);
	Vector3D j(a2.x, a2.y, a2.z);
	Vector3D k(a3.x, a3.y, a3.z);
	Vector3D l(a4.x, a4.y, a4.z);
	return Dihedral(i, j, k, l);
}


#endif
