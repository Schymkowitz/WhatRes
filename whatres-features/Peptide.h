// SequenceAnalysis, Copyright Switch Lab 2017
// Joost Schymkowitz

#ifndef PEPTIDE_H
#define PEPTIDE_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "Residue.h"
#include "vector3D.h"

using namespace std;


class Peptide
{
public:
	vector<Residue> reslist;
	map<int, int> res_contacts;
	vector<int> contact_dists;
	Vector3D center_of_gravity;
	
	char chain;
	map<string, char> three_to_one;
	double ContactOrder;
	double FoldDelay_Max;			// FoldDelay based on the most c-terminal interaction 
	double FoldDelay_Median;		// FoldDelay based on the median interactor
	double FoldDelay_Mean;			// FoldDelay based on the mean interactor
	double FoldDelay_80;			// FoldDelay based on 80% of the interactions formed
	
	Peptide(char x):chain(x), ContactOrder(0.), FoldDelay_Max(0.), FoldDelay_Median(0.), FoldDelay_Mean(0.), FoldDelay_80(0.)
		 { 
 	 		three_to_one.insert(pair<string,char>("ALA",'A'));
 	 		three_to_one.insert(pair<string,char>("GLY",'G'));
 	 		three_to_one.insert(pair<string,char>("LEU",'L'));
 	 		three_to_one.insert(pair<string,char>("VAL",'V'));
 	 		three_to_one.insert(pair<string,char>("ILE",'I'));
 	 		three_to_one.insert(pair<string,char>("PRO",'P'));
 	 		three_to_one.insert(pair<string,char>("ARG",'R'));
 	 		three_to_one.insert(pair<string,char>("THR",'T'));
 	 		three_to_one.insert(pair<string,char>("SER",'S'));
 	 		three_to_one.insert(pair<string,char>("CYS",'C'));
 	 		three_to_one.insert(pair<string,char>("MET",'M'));
 	 		three_to_one.insert(pair<string,char>("LYS",'K'));
 	 		three_to_one.insert(pair<string,char>("GLU",'E'));
 	 		three_to_one.insert(pair<string,char>("GLN",'Q'));
 	 		three_to_one.insert(pair<string,char>("ASN",'N'));
 	 		three_to_one.insert(pair<string,char>("ASP",'D'));
 	 		three_to_one.insert(pair<string,char>("TRP",'W'));
 	 		three_to_one.insert(pair<string,char>("TYR",'Y'));
 	 		three_to_one.insert(pair<string,char>("PHE",'F'));
 	 		three_to_one.insert(pair<string,char>("HIS",'H'));
		 }

	inline void Add_Residue(const Residue& res)
	{
		reslist.push_back(res);
	}
	
	inline void Set_Chain(char X)
	{
		chain=X;
	}
	
	inline int Rescount()
	{
		return reslist.size()-1;
	}
	
	inline void Add_Atom_to_lastrest(Atom& A)
	{
		reslist[reslist.size()-1].Add_Atom(A);
	}
	
	char convert_three(const string& three) const; 
	string GetSequence() const;
	void SetDihedrals();
	void SetCenterofGravity();
	void Translate(const double dx, const double dy, const double dz);
	void Rotate(const double theta, const Vector3D& unit_vector);
	void TruncateWindowToAlanine(int start, int len);

	void PrintDihedrals(ofstream& out, const string& pdb);
	void WriteContacts(ofstream& out, const string& pdb);
	void WritetoPDB(ofstream& o);
	void WritetoPDB(std::ostream& o);
};


inline ostream& operator<<(ostream& os, const Peptide& pep)  
{  
    os << pep.chain << endl;
	vector<Residue>::const_iterator it = pep.reslist.begin(), ite = pep.reslist.end();
	for(; it!=ite; it++) os << (*it);
    return os;  
} 

#endif
