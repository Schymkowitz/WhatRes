// SequenceAnalysis, Copyright Switch Lab 2020
// Joost Schymkowitz

#include "Peptide.h"

using namespace std;

Vector3D _rotate(const Vector3D& A, const double theta, const Vector3D& unit_vector)
{
	// Euler-Rodrigues' rotation formula
	// v_rot = v.cosØ + (k X v).sinØ + k(k.v)(1-cosØ)
	// where Ø is angle 
	// v is vector to rotate
	// k is unit vector around which to rotate, eg. (1,0,0)
	// X i cross product (in vector3D implemented as %)
	// . is dot product (in vector3D implemented as *)

	Vector3D v_rot, v_help1, v_help2, v_help3;
	
	v_help1 = cos(theta) * A;
	v_help2 = sin(theta) * (A % unit_vector);
    double helper = unit_vector * A;
    v_help3 = helper * unit_vector;
    v_help3 = (1 - cos(theta)) * v_help3;
	
	v_rot = v_help1 + v_help2 + v_help3;
	
	return v_rot;
	
}

void Peptide::TruncateWindowToAlanine(int start, int len)
{
    for(int i=0;i<len;i++)
    {
        int r = start + i;
        if(r < 0 || r >= (int)reslist.size()) continue;
        reslist[r].TruncateToAlanine();
    }
}

void Peptide::Rotate(const double theta, const Vector3D& unit_vector)
{
	vector<Residue>::iterator it = reslist.begin(), ite = reslist.end();
	for(; it != ite; it++)
	{
		vector<Atom>::iterator at = (*it).atoms.begin(), ate = (*it).atoms.end();
		for(; at != ate; at++)
		{
//			cout << *at;
            Vector3D v((*at).x, (*at).y, (*at).z);
			Vector3D v_rot = _rotate(v, theta, unit_vector);
			(*at).x = v_rot.x;
			(*at).y = v_rot.y;
			(*at).z = v_rot.z;
//			cout << *at << endl;
		}
	}	
}

void Peptide::Translate(const double dx, const double dy, const double dz)
{
	vector<Residue>::iterator it = reslist.begin(), ite = reslist.end();
	for(; it != ite; it++)
	{
		vector<Atom>::iterator at = (*it).atoms.begin(), ate = (*it).atoms.end();
		for(; at != ate; at++)
		{
	//		cout << *at;
			(*at).x += dx;
			(*at).y += dy;
			(*at).z += dz;
	//		cout << *at << endl;
		}
	}
}

void Peptide::SetCenterofGravity()
{
	double x(0.), y(0.), z(0.), beancount(0.);
	vector<Residue>::iterator it = reslist.begin(), ite = reslist.end();
	for(; it != ite; it++)
	{
		vector<Atom>::iterator at = (*it).atoms.begin(), ate = (*it).atoms.end();
		for(; at != ate; at++)
		{
			x += (*at).x;
			y += (*at).y;
			z += (*at).z;
			beancount++;
		}
	}
	x /= beancount;
	y /= beancount;
	z /= beancount;
	
	center_of_gravity = Vector3D(x, y, z);
}

char Peptide::convert_three(const string& three) const
{
	char ret = 'X';
      map<string, char>::const_iterator it = three_to_one.find(three);
	if(it == three_to_one.end()) return ret;
	else return (*it).second;
	return ret;
}
	
string Peptide::GetSequence() const
{
	stringstream helper;
	vector<Residue>::const_iterator it = reslist.begin(), ite = reslist.end();
	for(; it!=ite; it++) helper << convert_three((*it).residue_name);
      return helper.str();
}
	
void Peptide::SetDihedrals()
{
	// Every residue: phi and omega need the previous residue, psi the next one.
	// (The loop used to run 1..size-2, which left the first residue's psi and
	// the last residue's phi/omega at -999 although they are defined.)
	int size = reslist.size();
	for(int i=0; i < size; i++)
	{
		// phi
		if(i > 0 && reslist[i-1].checkC() && reslist[i].checkN() && reslist[i].checkCA() && reslist[i].checkC())
			reslist[i].phi = Dihedral(reslist[i-1].getC(), reslist[i].getN(), reslist[i].getCA(), reslist[i].getC());
		// psi
		if(i+1 < size && reslist[i].checkN() && reslist[i].checkCA() && reslist[i].checkC() && reslist[i+1].checkN())
			reslist[i].psi = Dihedral(reslist[i].getN(), reslist[i].getCA(), reslist[i].getC(), reslist[i+1].getN());
		// omega
		if(i > 0 && reslist[i-1].checkCA() && reslist[i-1].checkC() && reslist[i].checkN() && reslist[i].checkCA())
			reslist[i].omega = Dihedral(reslist[i-1].getCA(), reslist[i-1].getC(), reslist[i].getN(), reslist[i].getCA());
	}
}

void Peptide::PrintDihedrals(ofstream& out, const string& pdb)
{
	int size = reslist.size();
	for(int i=1; i < size-1; i++)
	{
		out << pdb << "\t" << reslist[i].residue_name << "\t" << reslist[i].chain << "\t" << i << "\t" << reslist[i].phi << "\t" << reslist[i].psi << "\t" << reslist[i].omega << endl;
	}
}

void Peptide::WriteContacts(ofstream& out, const string& pdb)
{
	out << "pdb\tresidue_name\tchain\tresnum\tmax_interactor\tContactOrder\n";
	int size = reslist.size();
	for(int i=0; i < size-1; i++)
	{
		out << pdb << "\t" << reslist[i].residue_name << "\t" << reslist[i].chain << "\t" << i << "\t" << reslist[i].max_interactor << "\t" << reslist[i].ContactOrder << endl;
	}
}


	
void Peptide::WritetoPDB(ofstream& o)
{
	int sizer = reslist.size();
	for(int i=0; i < sizer; i++)
	{
		reslist[i].WritetoPDB(o);
	}
}

void Peptide::WritetoPDB(std::ostream& o)
{
    int sizer = reslist.size();
    for(int i=0; i < sizer; i++)
    {
        reslist[i].WritetoPDB(o);
    }
}
