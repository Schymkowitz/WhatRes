// SequenceAnalysis, Copyright Switch Lab 2020
// Joost Schymkowitz

#include "PDB_dist_command.h"
#include "PDB_Reader.h"
#include "PDB_File.h"

#include <cmath>
#include <sstream>     
#include <ctime>


void PDB_dist_command::SetOptionsRequired()
{
	
}

void PDB_dist_command::Execute()
{
    cout << "TTTesting..." << endl;

    ifstream pdblist(Option.inputfiles[0]);
	ofstream outer("PSX_PDB_dist.txt"), outer2("PSX_PDB_dihedrals.txt");
	outer2 << "pdb\tres\tchain\tnumber\tphi\tpsi\tomega\n";
	string currentpdb;
    vector<string> pdbs;
    while(!pdblist.eof())
    {
    	string line;
     	GetLine(pdblist,line);
    	if(line.size()==0) continue;
		pdbs.push_back(line);
	}
	pdblist.close();
	
	int pdblistsize = pdbs.size();
	for(int i=0; i < pdblistsize; i++)
	{
		AnalysePDB(pdbs[i], outer2);
	}
	
	// write output
	HH.Print(outer);
	
	cout << "number of histograms: " << HH.Size() << endl;
	outer.close();
	outer2.close();
	cout << pdblistsize << "files read and processed." << endl;
	
}
	
	
void PDB_dist_command::AnalysePDB(const string& thispdb, ofstream& dihedout)
{	
	// extracts atomic distance distributions from structures
	
    cout << thispdb << endl;
    PDB_Reader reader;
    PDB_File pdb;
    reader.Readpdb(thispdb, pdb);

	pdb.PrintDihedrals(dihedout);

	int pepnum = pdb.peptides.size();
	vector< Residue > allres;
    for(int e=0; e < pepnum; e++)
    {
		Peptide& thispep = pdb.peptides[e];
		int resnum = thispep.reslist.size();
		for(int q=0; q < resnum; q++)
		{
			Residue& thisres = thispep.reslist[q];
			allres.push_back(thisres);
		}
	}
	
	int allresnum = allres.size();
	
	for(int i = 0; i < allresnum; i++)
	{
		Residue& firstres = allres[i];
		int firstatomnum = firstres.atoms.size();
		
		for(int e = i+1; e < allresnum; e++)
		{
			Residue& secondres = allres[e];
//			cout << firstres.chain << "\t" << secondres.chain << "\t" << firstres.number_in_pdb << "\t" << secondres.number_in_pdb << endl;
			if(secondres.chain == firstres.chain && abs(firstres.number_in_pdb - secondres.number_in_pdb) < 2) 
			{
//				cout << "abort" << "\t" << firstres.number_in_pdb - secondres.number_in_pdb << endl;
				continue;
			}
			int secondatomnum = secondres.atoms.size();
			
			for(int o = 0; o < firstatomnum; o++)
			{
				Atom& a1 = firstres.atoms[o];
				for(int q = 0; q < secondatomnum; q++)
				{
					Atom& a2 = secondres.atoms[q];
					double _dist = pdb.dist(a1, a2);
					if(_dist < DIST_CUT)
					{ 
						stringstream handle1, handle2;
						handle1 << firstres.residue_name << "." << a1.atom_name;
						string type1 = a1.GetAtomType(handle1.str());
						if(type1 == "nonsense") continue;
						handle2 << secondres.residue_name << "." << a2.atom_name;
						string type2 = a2.GetAtomType(handle2.str());
						if(type2 == "nonsense") continue;
						
						HH.Add(type1, type2, _dist);
					}
				}
			}
		}
	}
	
}
