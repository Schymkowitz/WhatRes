// SequenceAnalysis, Copyright Switch Lab 2017
// Joost Schymkowitz

#ifndef ATOM_H
#define ATOM_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;


class Atom
{
private:
	static std::map<std::string, std::string>& AtomTypes();
	static void EnsureAtomTypesLoaded(const std::string& style);

public:
	string atom_name;
	double x;
	double y;
	double z;
	double bfactor;
	
    Atom(const string& AN, double X, double Y, double Z, double BF):atom_name(AN), x(X), y(Y), z(Z), bfactor(BF)
	  {
	  }
	  
   static void LoadAtomTypes(const string& style);
   string GetAtomType(const string& handle);
   void PutPDBline(ostream &o, double BF, int res_index, int atom_index, char chain, string &AAtype);
   
};


inline ostream& operator<<(ostream& os, const Atom& at)  
{  
    os << at.atom_name << "\t" << at.x << "\t" << at.y << "\t" << at.z << "\t" << at.bfactor << endl;
    return os;  
} 

#endif
