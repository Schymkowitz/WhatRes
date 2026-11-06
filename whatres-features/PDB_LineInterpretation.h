// SequenceAnalysis, Copyright Switch Lab 2017
// Joost Schymkowitz


#ifndef PDBLINEINTERPRETATION_H
#define PDBLINEINTERPRETATION_H

#include <string>
#include <iostream>
#include <cstring>
#include <cstdlib>
using namespace std;

class PDB_LineInterpretation
{
public:
	string keyword;
	string atom;
	string amino;
	string abc;
	char chain;
	int index;
	int resnumber;
	string pdbseqlet;
	string resnumbackup; // read residue numbers as strings for checks
	double x;
	double y;
	double z;
	double b_factor;
	bool proton;
	bool water;
	bool valid;

	void output_line()
	{
        cout << keyword << " " << atom << " " << amino << " " << abc << " " << index << " " << resnumber << " " << x << " " << y << " " << z << " " << endl;
    }

	PDB_LineInterpretation():index(0), resnumber(0), x(0.0), y(0.0), z(0.0), b_factor(0.0), proton(false), water(false), valid(false) { }

	void clean(std::string &s)
	{
		size_t i = 0;
		i = s.find(" ");
		while (i!=string::npos)
		{
			s = s.erase(i, 1);
			i = s.find(" ");
		}
	}

	std::string GetS(std::string s1, int start, int num)
	{
		string s2;
		s2.assign(s1, start, num);
		clean(s2);
		return s2;
	}

	int GetI(std::string s1, int start, int num)
	{
		std::string s2;
		int i = 00;
		s2 = GetS(s1, start, num);
		i = atoi(s2.c_str());
		return i;
	}

	double GetF(string s1, int start, int num)
	{
		string s2;
		double f = 0.0;
		s2 = GetS(s1, start, num);
		f = atof(s2.c_str());
		return f;
	}

	bool standardise_atom();

	void check_atom()
	{
		string s;
		string c;
		water = proton = false;
		if (amino == "HOH")
		{
			water = true;
			return;
		}
		else
		{
			// Hydrogen, or deuterium (neutron structures: D, DG1, DD21, 1DB...);
			// no standard heavy-atom name starts with D.
			s.assign(atom, 0, 1);
			if (s == "H" || s == "D") { proton = true; return; }
			c.assign(atom, 1, 1);
			if ((s == "1" || s == "2" || s == "3" || s == "4") && (c == "H" || c == "D"))
				proton = true;
		}
	}

	void Atom_reader(string &line)
	{
		atom = GetS(line, 12, 4);
		if (atom[0]=='H' || (atom.size()>2 && (atom[1]=='H' && atom[2]=='H'))) { keyword="PROTON"; return; }
		amino = GetS(line, 17, 3);
		abc = GetS(line, 16, 1);
		chain = ' ';
		string _chain = GetS(line, 21, 1);
		if (_chain.size() >0) chain = _chain[0];
		resnumber = GetI(line, 22, 4);
		pdbseqlet = GetS(line, 26, 1);
		resnumbackup = GetS(line, 22, 5);
		index = GetI(line, 6, 5);
		x = GetF(line, 30, 8);
		y = GetF(line, 38, 8);
		z = GetF(line, 46, 8);
		if (line.size() >=66 )
			b_factor = GetF(line, 60, 6);
		else
			b_factor = 0.0;
		check_atom();
		valid = standardise_atom();
	}

	void HetAtom_reader(string &line)
	{
		atom = GetS(line, 12, 4);
		amino = GetS(line, 17, 3);
		abc = GetS(line, 16, 1);
		chain = ' ';
		string _chain = GetS(line, 21, 1);
		if (_chain.size() >0) chain = _chain[0];
		resnumber = GetI(line, 22, 4);
		pdbseqlet = GetS(line, 26, 1);
		resnumbackup = GetS(line, 22, 5);
		index = GetI(line, 6, 5);
		x = GetF(line, 30, 8);
		y = GetF(line, 38, 8);
		z = GetF(line, 46, 8);
		if (line.size() >=66 )
			b_factor = GetF(line, 60, 6);
		else
			b_factor = 0.0;
		check_atom();
		if (!water && !proton)
			valid = standardise_atom();
	}


	void Remark_reader(string line)
	{
		index = GetI(line, 7, 4);
		if (index == 350 && line.size() > 20)
		{
			atom = GetS(line, 11, 9);
			atom.assign(atom, 0, 5);
			if (atom == "BIOMT")
			{
				resnumber = GetI(line, 21, 3);
				x = GetF(line, 24, 10);
				y = GetF(line, 34, 10);
				z = GetF(line, 44, 10);
			}
		}
	}

	void interpret(string &line)
	{
		std::string temp;
		valid = true;
		keyword.assign(line, 0, 6);
		clean(keyword);
		if (keyword == "ATOM")
			Atom_reader(line);
		else if (keyword == "HETATM")
		{
			HetAtom_reader(line);
			if (amino=="TPO" ||amino=="SEP" ||amino=="PTR")
			{
				keyword="ATOM";
				Atom_reader(line);
			}
		}
		else if (keyword == "REMARK")
			Remark_reader(line);
	}
};

#endif

