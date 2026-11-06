// SequenceAnalysis, Copyright Switch Lab 2019
// Joost Schymkowitz

#include "PDB_File.h"
#include <math.h>
#include <algorithm>
#include <numeric>


bool PDB_File::CheckChain(char c)
{
	vector<char>::iterator it = find(chain_names.begin(), chain_names.end(), c);
	return (it != chain_names.end());
}

void PDB_File::WriteChainPDBFormat(ofstream& out, char c)
{
	out << "REMARK SwitchLab PDB File (C) 2020\n";
	
	int sizer = peptides.size();
	for(int i=0; i < sizer; i++)
	{
		if(chain_names[i] == c)
		{
			peptides[i].WritetoPDB(out);
			out << "TER\n";
		}
	}
	out << "ENDMDL\n";	
}



void PDB_File::WritePDBFormat(ofstream& out)
{
	out << "REMARK SwitchLab PDB File (C) 2020\n";

	// Pipeline provenance accumulated by ingest stages.
	for (const auto& ln : provenanceLines)
		out << "REMARK 999 " << ln << "\n";

	// Re-emit experimental metadata so the inspection PDB carries it.
	if (!experimentalMethod.empty())
		out << "EXPDTA    " << experimentalMethod << "\n";
	if (resolutionA >= 0.0) {
		char buf[64];
		snprintf(buf, sizeof(buf), "REMARK   2 RESOLUTION. %6.2f ANGSTROMS.\n", resolutionA);
		out << buf;
	}

	int sizer = peptides.size();
	for(int i=0; i < sizer; i++)
	{
		peptides[i].WritetoPDB(out);
		out << "TER\n";
	}
	out << "ENDMDL\n";

}

void PDB_File::PrintDihedrals(ofstream& out)
{
	int sizer = peptides.size();
	for(int i=0; i < sizer; i++)
	{
		peptides[i].PrintDihedrals(out, pdbName);
	}
}

double PDB_File::dist(Atom& a1, Atom& a2)
{
	double x = pow(a1.x - a2.x, 2);
	double y = pow(a1.y - a2.y, 2);
	double z = pow(a1.z - a2.z, 2);
	
	double _dist = sqrt(x + y + z);
	return _dist;
}

void PDB_File::SetContactOrder(double dist_cutoff, double seqdist_cutoff)
{
	bool _DEBUG = false;
	int sizer = peptides.size();
	for(int e = 0; e < sizer; e++)
	{
		Peptide& P = peptides[e];
		int sozer = P.reslist.size();

		for(int i = 0; i < sozer; i++)
		{
			Residue& Resi = P.reslist[i];
			if(Resi.residue_name == "GAP") continue;
			int atcounti = Resi.atoms.size();

			for(int j=0; j < sozer; j++)
			{
				double seqdist = abs(j-i);
				if( seqdist <  seqdist_cutoff) continue;

				Residue& Resj = P.reslist[j];
				if(Resj.residue_name == "GAP") continue;
				int atcountj = Resj.atoms.size();
				bool done = false;
				
				for(int ai = 0; ai < atcounti && !done; ai++)
				{
					Atom& Ai = Resi.atoms[ai];
                    if(Ai.atom_name[0] == 'H') continue; // skip hydrogens

					for(int aj = 0; aj < atcountj && !done; aj++)
					{
						Atom& Aj = Resj.atoms[aj];
                        if(Aj.atom_name[0] == 'H') continue; // skip hydrogens

						double aadist = dist(Ai, Aj);
						if(aadist > 20) done = true;  // for speed, if any interatomic distance is > 20 then we assume these are no other atom contacts for these residues
						
						if(aadist < dist_cutoff)
						{
							if(_DEBUG)
							{
								cout << Resi.residue_name << Resi.number_in_pdb << " " << Ai.atom_name << "  -  "
                                << Resj.residue_name << Resj.number_in_pdb << " " << Aj.atom_name << "\n";
							}
							// atom contact detected between Resi and Resj
							Resi.seqdist.push_back(seqdist);
							Resj.seqdist.push_back(seqdist);
							Resi.abs_int.push_back(j);
							Resj.abs_int.push_back(i);
							if(Resi.max_interactor < j) Resi.max_interactor = j;
							if(Resj.max_interactor < i) Resj.max_interactor = i;
							if(i < j)
							{
								P.contact_dists.push_back(j-i);
								P.res_contacts.insert(pair<int,int>(i,j)); 
							}
						}
					}
				}
			}
		}
	//	cout << "res\tsize\tmax\tmean\tmedian\tpos80%\tdel80\n";
		for(int i = 0; i < sozer; i++)
		{
            Residue& Resi = P.reslist[i];
			if(Resi.residue_name == "GAP") continue;
            vector<int>& v = Resi.seqdist;
			double seqdistsum(0.);
			if(v.size() > 0)
			{
				seqdistsum = accumulate( v.begin(), v.end(), 0.0);
				Resi.ContactOrder = seqdistsum / v.size() ;
			}
			else
			{
				Resi.ContactOrder = 0.;
			}
			if(Resi.max_interactor > i)
			{
				// max interactor is c-terminal, count the delay
				P.FoldDelay_Max += (Resi.max_interactor - i);
			}
			size_t intcount = Resi.abs_int.size();
			if(intcount < 1) continue;
			
			sort(Resi.abs_int.begin(), Resi.abs_int.end());
			Resi.interactor_mean = floor( accumulate(Resi.abs_int.begin(), Resi.abs_int.end(), 0.) / Resi.abs_int.size() ) - i;
            int median_pos = floor(Resi.abs_int.size()/2);
            Resi.interactor_median = Resi.abs_int[median_pos] - i;
			int pos80 = floor(0.8 * Resi.abs_int.size());
            Resi.interactor_80 = Resi.abs_int[pos80] - i;
	
			if(Resi.interactor_mean > 0) P.FoldDelay_Median += Resi.interactor_mean;
			if(Resi.interactor_median > 0) P.FoldDelay_Mean += Resi.interactor_median;
			if(Resi.interactor_80 > 0) P.FoldDelay_80 += Resi.interactor_80;
            
			if(Resi.interactor_median > 500)
			{
				cout << i << "\t" << Resi.abs_int.size() << "\t" << Resi.max_interactor << "\t" << Resi.interactor_mean << "\t" << Resi.interactor_median << "\t" << Resi.interactor_80 << endl;
			}
		}
		P.FoldDelay_Max /= sozer;
		P.FoldDelay_Median /= sozer;
		P.FoldDelay_Mean /= sozer;
		P.FoldDelay_80 /= sozer;		
		
		double peptide_seqdistsum(0.), peptide_numcontacts(P.contact_dists.size());
		vector<int>::iterator it = P.contact_dists.begin(), ite = P.contact_dists.end();
		for(; it != ite; it++)
		{
			peptide_seqdistsum += (*it);
		}
		P.ContactOrder = peptide_seqdistsum / sozer / peptide_numcontacts;
		cout << "\t" << P.chain << "\t" << P.ContactOrder << "\t" << peptide_seqdistsum << "\t" << peptide_numcontacts << "\t" << sozer << endl;
	}
}

void PDB_File::Sequence_on_screen() 
{
	for(int i=0; i<peptides.size(); i++)
	{
		cout << pdbName << "_" << peptides[i].chain << endl;
		cout << sequences[i] << endl;
		cout << endl;
	}
}

bool PDB_File::Set_currentP(char a)
{
	for(int i=0; i < peptides.size(); i++)
	{
		if(peptides[i].chain == a)
		{
			Set_currentP(i);
			return true;
		}
	}
	return false;
}

int PDB_File::identify_chain(const char c)
{
	int pos = -1;
	for( int i=0; i < chain_names.size(); i++)
	{
		if(chain_names[i] == c)
		{
			pos = i;
			break;
		}
	}
	if(pos==-1)
	{
		cout << "chain identifier does not exist in pdb file " << endl;
	}
	return pos;		
}
int PDB_File::Get_Pos_from_PDBnumber(int pdbnum, const char c)
{
	int ch = identify_chain(c);
	if(ch==-1)
		return -1;
	Peptide& thispep = peptides[ch];
	for(int i=0; i< thispep.reslist.size(); i++)
	{
		Residue& thisres = thispep.reslist[i];
		if(thisres.number_in_pdb == pdbnum)
		{
			return i;
		}
	}
	cout << "specified residuenumber not found in this pdb" << endl;
	return -1;
}

string PDB_File::GetChainSeq(const char c)
{
	string helper;
	int pos = identify_chain(c);
	if(pos != -1)
		return sequences[pos];
	return helper;	
}

string PDB_File::GetChainSeq(int pepnum)
{
	string hlpr;
	if(pepnum < 0 || pepnum > sequences.size() -1 )
	{
		cout << "impossible peptide number in pdb requested: " << pepnum << endl;
		return hlpr;
	}
	return sequences[pepnum];
}

