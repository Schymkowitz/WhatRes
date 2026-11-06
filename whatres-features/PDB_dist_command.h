// SequenceAnalysis, Copyright Switch Lab 2010
// Joost Schymkowitz

#ifndef PDB_DIST_COMMAND_H
#define PDB_DIST_COMMAND_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string.h>
#include <map>

#include "commands.h"
//#include "translation.h"
//#include "blast.h"
#include "number_histogram.h"

using namespace std;

class histohelper
{
	const double MINDIST = 0.0;
	const double MAXDIST = 9.0;
	const int BINNUM = 60;
	map<string, NumberHistogram> histos;
 	
public:
	histohelper()
	{}
		
	void Add(const string& part1, const string& part2, double dist)
	{
		// ensure the atom names are always in the same order
		stringstream handle;
        if(part1 < part2) handle << part1 << "_" << part2;
		else handle << part2 << "_" << part1;

		map<string, NumberHistogram>::iterator it = histos.find(handle.str());

		if(it == histos.end())
		{
			// no hist found for this handle
			// create new histo and stick in list
			NumberHistogram NH(BINNUM, MINDIST, MAXDIST);
			NH.Add(dist);
			histos.insert(pair<string, NumberHistogram>(handle.str(), NH));
		}
		else
		{
			// already started a hist for this handle
			(*it).second.Add(dist);
		}
	}
	
	int Size()
	{
		return histos.size();
	}
	
	void Print(ofstream& outer)
	{
		map<string, NumberHistogram>::iterator it = histos.begin(), ite = histos.end();
		double binwidth = (MAXDIST - MINDIST) / BINNUM;
		outer << "dist\t";
		for(double i = MINDIST; i < MAXDIST; i += binwidth)
			outer << i + 0.5 * binwidth << "\t";
		outer << "\n";
		for(; it != ite; it++)
		{
			(*it).second.Totalise();
			(*it).second.Normalise();
			outer << (*it).first << "\t";
			(*it).second.PrintFractions(outer);
			outer << "\t" << (*it).second.Get_totalcount() << "\n";
		}
	}
};

class PDB_dist_command: public sequence_command
{
	const double DIST_CUT = 9.0;
	histohelper HH;
	 
public:
	PDB_dist_command(Option_class& O):sequence_command("PDB_dist", "Analyse pairwise atom distances in PDBs", O)
	{
	}

	void Execute() override;
	void SetOptionsRequired() override;
	void AnalysePDB(const string& thispdb, ofstream& dihedout);

	std::string HelpText() override
	{
		return R"(Build reference distributions of pairwise atom distances and backbone
dihedrals from a set of PDB structures.

What it does
  Reads a list of PDB files. For each it records, per residue, the
  phi/psi/omega backbone dihedrals, and accumulates atom-pair distance
  histograms: every non-adjacent residue pair contributes each atom pair
  closer than the 9 A cutoff, binned by atom-type pair.

  This is a structural-statistics tool. It overlaps with structure_survey
  but measures different things - the pairwise atom-type distance
  distributions collected here are not produced by structure_survey. It
  provided the pairwise atom-distance analysis in van der Kant et al.,
  Structure (2022).

Input
  -in   a .txt file listing one PDB path per line.

Output
  PSX_PDB_dist.txt        atom-type-pair distance histograms (60 bins,
                          0-9 A), normalised, with total counts
  PSX_PDB_dihedrals.txt   per-residue phi / psi / omega)";
	}
};

#endif
