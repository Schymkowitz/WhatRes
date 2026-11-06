// SequenceAnalysis, Copyright Switch Lab 2017
// Joost Schymkowitz

#include "PDB_Reader.h"
#include "LineReader.h"
int PDB_Reader::Readpdb(const string &naam, PDB_File& pdb)
{
    int model_count=0;
	const int starter = -999;
	char helpchain('$');
	std::string gap("GAP"), sum="";
	int Rcount = starter, Pcount = 0;
	string line;
	int helpres = starter;
	string helpresnumbackup = "XXXTTXXX";
	// altLoc collapse state: the conformer label kept for the current residue
	string altResKey = "", altKeep = "";

	bool PDBok = true;
	std::string unknown_hetatm("Unknown heteroatoms in pdb: ");
	vector<string> unknown_list;
	int hlp;
	pdb.pdbName = naam;

	// Quick MODEL-count pass.  PDB_Reader's main loop stops at the first
	// ENDMDL (taking only model 1 for NMR ensembles), so MODEL records
	// beyond the first wouldn't otherwise be seen.  This count exposes
	// the ensemble signal to downstream tools via PDB_File::nModels.
	{
		std::ifstream pre(naam.c_str());
		if (pre.good()) {
			std::string preLine;
			while (std::getline(pre, preLine)) {
				if (!preLine.empty() && preLine.back() == '\r') preLine.pop_back();
				if (preLine.compare(0, 6, "MODEL ") == 0) ++pdb.nModels;
			}
		}
	}

	if(OpenPdb(naam) < 0) return -1;

	int mdlcount=0;
	if(model_count!=0) // for nmr
	{
		while (!pdbfile.eof())
		{
			GetLine(pdbfile, line);
			LIP.interpret(line);
			if(LIP.keyword=="ENDMDL") 
			{
				mdlcount++;	  
				if(mdlcount==model_count) break;
			}
		}
	}
	
	if(pdbfile.eof()) return -1;
	while (!pdbfile.eof())
	{
		GetLine(pdbfile, line);

		// Capture header records that LIP would otherwise skip / mark
		// invalid.  Downstream tools (BiomtExpander, CrystalLayerExpander)
		// read these from PDB_File rather than re-opening the source file.
		if (line.compare(0, 6, "CRYST1") == 0 && pdb.cryst1Line.empty())
			pdb.cryst1Line = line;
		else if (line.compare(0, 11, "REMARK 290 ") == 0)
			pdb.remark290Lines.push_back(line);
		else if (line.compare(0, 11, "REMARK 350 ") == 0)
			pdb.remark350Lines.push_back(line);
		else if (line.compare(0, 6, "EXPDTA") == 0 && pdb.experimentalMethod.empty()) {
			// EXPDTA contents start in column 11 (1-based); trim trailing space.
			std::string m = line.size() > 10 ? line.substr(10) : std::string();
			while (!m.empty() && (m.back() == ' ' || m.back() == '\t' || m.back() == '\r'))
				m.pop_back();
			while (!m.empty() && m.front() == ' ') m.erase(m.begin());
			pdb.experimentalMethod = m;
		}
		else if (line.compare(0, 10, "REMARK   2") == 0 && pdb.resolutionA < 0.0) {
			// "REMARK   2 RESOLUTION.    X.XX ANGSTROMS.": extract the number.
			std::size_t p = line.find("RESOLUTION");
			if (p != std::string::npos) {
				std::size_t i = p + 10;
				while (i < line.size() && (line[i] == ' ' || line[i] == '.' || line[i] == ':')) ++i;
				if (i < line.size() && (std::isdigit((unsigned char)line[i]) || line[i] == '.')) {
					pdb.resolutionA = std::atof(line.c_str() + i);
				}
			}
		}

		LIP.interpret(line);

		if (LIP.keyword=="ANISOU")
			continue;
		if(LIP.keyword=="ENDMDL")
		{
			model_count++;
			break;
		}
		if (LIP.proton) continue;
		else if (!LIP.valid)
		{
			unknown<<line<<endl;
			continue;
		}

		// Alternate-conformer (altLoc) collapse.  The altLoc indicator
		// (column 17) is parsed into LIP.abc.  A structure that deposits
		// multiple side-chain rotamers tags them A/B/C...; loading every
		// one inflates the atom count and creates spurious intra-residue
		// clashes that trip the geometry gate (e.g. 5AEF, 9KTZ).  Keep
		// only one conformer per residue: atoms with a blank altLoc plus
		// the FIRST altLoc label seen in that residue (normally "A").
		// Keeping literally "A" lost every atom of residues whose
		// conformers are labelled only B/C/D (e.g. 9HQ6 Arg67: C and D).
		if (LIP.keyword == "ATOM" || LIP.keyword == "HETATM")
		{
			// residue = chain + columns 23-27 (number + insertion code); not the
			// residue name, so microheterogeneity (A=ALA, B=SER) still keeps one
			string resKey = string(1, LIP.chain) + "|" + LIP.resnumbackup;
			if (resKey != altResKey) { altResKey = resKey; altKeep = ""; }
			if (!LIP.abc.empty())
			{
				if (altKeep.empty()) altKeep = LIP.abc;
				if (LIP.abc != altKeep) continue;
			}
		}

		while (LIP.keyword == "REMARK" && LIP.index == 350)
		{
			GetLine(pdbfile, line);
			// Capture continuation REMARK 350 lines that we advance through.
			if (line.compare(0, 11, "REMARK 350 ") == 0)
				pdb.remark350Lines.push_back(line);
			LIP.interpret(line);
		}
		if ((LIP.keyword == "ter" || LIP.keyword == "TER"))
		{
// 			cout << LIP.keyword<<"\t"<<LIP.chain<<"\t"<<helpchain<<endl;
// 			helpchain = '$';
			continue;
		}
		if (LIP.keyword == "HETATM" && (LIP.amino=="PTR" || LIP.amino=="TPO" || LIP.amino=="SEP"
			|| LIP.amino=="MLY" || LIP.amino=="MLZ" || LIP.amino=="M3L" || LIP.amino=="6MA" || LIP.amino=="5CM" || LIP.amino=="HYP" || LIP.amino=="CSE"))
		{LIP.keyword = "ATOM" ;  }
		if (LIP.keyword == "ATOM" && (LIP.amino=="ZN" || LIP.amino=="MG" || LIP.amino=="CA" || LIP.amino=="MN" || LIP.amino=="CU" || LIP.amino=="FE"))
		{
			LIP.keyword = "HETATM"; 
		}
		if (LIP.keyword == "ATOM" && LIP.atom=="O" && LIP.amino=="HOH")
		{
			LIP.keyword = "HETATM"; 
		}
		if (LIP.keyword == "HETATM")
		{
            // water, metal, etc -> ignored for now
		}

		if (LIP.keyword == "ATOM" && !LIP.proton && (sum != (LIP.amino + LIP.atom)))
		{
			if (PDBok)
				PDBok = LIP.valid;
			
			// if the chain is interrumpted by a chainless residue, then assign the chain name of the previous residue
			if(LIP.chain == ' ' && helpchain != '$') {
				LIP.chain = helpchain;
			}

			if (helpchain != LIP.chain)
			{
				MolNum++;
                Peptide P(LIP.chain);
                pdb.Add_Peptide(P);
				pdb.Set_lastP_as_current();
				Rcount = starter;
				helpres = starter;
			}
			Pcount = pdb.Pepcount();
			if (helpres != LIP.resnumber || helpresnumbackup != LIP.resnumbackup)
			{
				if (Rcount == starter || LIP.resnumber - helpres <= 1)
				{
					hlp = Rcount+1;
					if (Rcount == starter)
						hlp = 0;
                    Residue R(LIP.keyword, LIP.amino, LIP.chain, LIP.resnumber);
					//if(LIP.pdbseqlet != "") thisRes->setPdbSeqLet(LIP.pdbseqlet);
                    pdb.Add_Res_to_currentP(R);
					Rcount = pdb.CurrentP_ResCount();
					if (helpchain != LIP.chain)
					{
                        /*
						if (thisRes->r_type == DNA)
							thisPep->SeMoleculeIdentity(M_DNA);
						if (thisRes->r_type == PROTEIN)
							thisPep->SeMoleculeIdentity(M_PROTEIN);
						if (thisRes->r_type == SMALLMOLECULE)
							thisPep->SeMoleculeIdentity(M_SMALLMOLECULE);
						if (thisRes->r_type == RNA)
							thisPep->SeMoleculeIdentity(M_RNA);
						if (thisRes->r_type == GAP)
							thisPep->SeMoleculeIdentity(M_GAP);
						if (thisRes->r_type == UNDEFINED)
							thisPep->SeMoleculeIdentity(M_UNDEFINED);
                         */
					}
				}
				else
				{
					for (int i = helpres+1; i < LIP.resnumber; i++) // add gap residues
					{
                        Residue R(LIP.keyword, gap, LIP.chain, i);
                        pdb.Add_Res_to_currentP(R);
						Rcount = pdb.CurrentP_ResCount();
					}
					
                    Residue R(LIP.keyword, LIP.amino, LIP.chain, LIP.resnumber);
                    pdb.Add_Res_to_currentP(R);
                    Rcount = pdb.CurrentP_ResCount();

					if (helpchain != LIP.chain)
					{
                        /*
						if (thisRes->r_type == DNA)
							thisPep->SeMoleculeIdentity(M_DNA);
						if (thisRes->r_type == PROTEIN)
							thisPep->SeMoleculeIdentity(M_PROTEIN);
						if (thisRes->r_type == SMALLMOLECULE)
							thisPep->SeMoleculeIdentity(M_SMALLMOLECULE);
						if (thisRes->r_type == RNA)
							thisPep->SeMoleculeIdentity(M_RNA);
						if (thisRes->r_type == GAP)
							thisPep->SeMoleculeIdentity(M_GAP);
						if (thisRes->r_type == UNDEFINED)
							thisPep->SeMoleculeIdentity(M_UNDEFINED);
                         */
					}
				}
			}
			if (Rcount == 0)
			{
                /*
				if (Options.Put_Nterm(thisRes->chain()))
					thisRes->setIsNterm(true);
				thisRes->setIsFirst(true);
                 */
			}
// 			if (( thisRes->ThreeLetter != "PRO" && thisRes->ThreeLetter != "HYP" ) && thisRes->isNterm() && LIP.atom == "N")
// 				LIP.atom = "tN";
            Atom A(LIP.atom, LIP.x, LIP.y, LIP.z, LIP.b_factor);
            pdb.Add_Atom_to_lastrest(A);
			helpchain = LIP.chain;
			helpres = LIP.resnumber;
			helpresnumbackup = LIP.resnumbackup;
			sum = LIP.amino + LIP.atom;
		}
	}

	for(int i=0; i<pdb.peptides.size(); i++)
	{
		pdb.sequences.push_back(pdb.peptides[i].GetSequence());
		pdb.chain_names.push_back(pdb.peptides[i].chain);
	}
    /*
	if(RA.peptide.size()==0) return -1;
	for (int i = 0; i < RA.peptide.size(); i++)
	{
		Residue *R=RA.peptide[i]->residue[RA.peptide[i]->residue.size()-1];
		if (Options.Put_Cterm(R->chain()))
			R->setIsCterm(true);
		R->setIsLast(true);
	}
	if (unknown_list.size() > 1000)
	{
		cerr << unknown_hetatm << endl;
		PDBok = false;
	}
    */
	ClosePdb();
	
	int pepsiz = pdb.peptides.size();
	for(int i =0; i < pepsiz; i++)
	{
		pdb.peptides[i].SetDihedrals();
	}

	if(!PDBok) return -1;
	
	
	return model_count;

}

