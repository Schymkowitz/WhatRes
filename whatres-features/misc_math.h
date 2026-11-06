// SequenceAnalysis, Copyright Switch Lab 2020
// Joost Schymkowitz
// This library is for simple math- and allocation- operations. It includes numrec ran2 function and as well as a c++ version of it 

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <cstring>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <iterator>

const double PI = 4.0*atan(1.0);
const double LOGPI=log(PI);
const double LOG2=log(2.0);
const double SMALL = 1.0e-15;
const double MINF = -1.0e+15;
const double INF = 1.0e+15;

#ifdef _MSC_VER
#pragma warning (disable : 4127)
#endif

using namespace std;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// permutation of a vector, stores the permuted vectors in a result vector of vectors
// takes as input a pointer to a vector, start position (starts with 0), the size of
// the vector and a pointer to a vector of pointers to vectors ;-)

template<typename T>
void permuteSize(vector<T>* v, const int start, const int n, vector<vector<T>*>* results)
{
	if (start == n-1)
	{
		vector<T>* newvec(new vector<T>);
		for (typename vector<T>::iterator k=v->begin(); k!=v->end(); k++)
		{
			T tmpvalue = *k;
			//cout << tmpvalue << ",";
			newvec->push_back(tmpvalue);
		}
		results->push_back(newvec);
		//cout << "test hier :"  << result->size() << endl;
		//cout << "]" << endl;
	}
	else
	{
		for (int i = start; i < n; i++)
		{
			//swap v[i] and v[start]
			T tmp1 = (*v) [i];
			T tmp2 = (*v) [start];
			(*v) [i] = tmp2;
			(*v) [start] = tmp1;
			permuteSize(v, start+1, n, results);
			//swap back again
			(*v) [start] = tmp2;
			(*v) [i] = tmp1;
		}
	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// permutation of all subsets of a given size (n), stores these in a vector of vectors
// uses the permuteSize function

template<typename T>
void permutationSize(vector<T>* all_values, vector<T> *selected_values, int n, vector<vector<T>*>* results)
{
	if (all_values->empty())
	{
		if (selected_values->size() == n)
		{
			// No more remaining values to be treated, print the permutation of this set
			//Loop until all permutations are generated.

			permuteSize(selected_values, 0, n, results);
		}
	}

	else
	{
		// Take the first element of all_values
		typename vector<T>::iterator it = all_values->begin();
		const T tmp = *it;

		all_values->erase(all_values->begin());

		// Call recursively without the first element added
		permutationSize(all_values, selected_values, n, results);
		// Call recursively with the first element added
		selected_values->insert(selected_values->end(), tmp);
		permutationSize(all_values, selected_values, n, results);
		// restore the sets
		it = selected_values->end(); it--;
		selected_values->erase(it);
		all_values->insert(all_values->begin(), tmp);
	}
}

//-------------
void newt_hess(double x[], int n, int *check, void (*vecfunc) (int, double [], double []), void (*Hessian) (int n, double x[], double fvec[], double **df), int NEWT_DEBUG=0);


inline void Addstrings(char *result, const char *s1, const char *s2)
{
	char *temp=(char *) malloc((strlen(s1)+strlen(s2)+1)*(sizeof(char)));
	strcpy(temp, s1);
	strcpy(temp+strlen(s1), s2);
	strcpy(result, temp);
	free(temp);
}


template<class T>
inline T Max(T *set, int nelem, int *where)
{
	int i;
	*where=0;
	T m=set[0];
	for (i=1; i<nelem; i++)
		if (set[i]>m) { m=set[i]; *where=i; }
	return m;
}

template<class T>
inline T Min(T *set, int nelem, int *where)
{
	int i;
	T m=set[0];
	*where=0;
	for (i=1; i<nelem; i++)
		if (set[i]<m) { m=set[i]; *where=i; }
	return m;
}

template<class T>
inline int InSet(T *set, T quan, int nset, int *elem)
{
	int i;
	for (i=0; i<nset; i++)
	{
		if (fabs(static_cast<double>(quan-set[i]))<SMALL) { *elem=i; return 1; }
	}
	*elem=-1;
	return 0;
}

template<class T>
inline void RemoveVal(T *set, T val, int fromset, int toset, int *ind, int *nind)
{
	int i=0, j=0;
	for (i=fromset; i<toset; i++)
	{
		if (fabs(set[i]-val)>SMALL) { ind[j]=i; j++; }
	}
	*nind=j;
}

template<class T>
inline  void PickSmallerThan(int *pickedindicies, int *npval, T *function, int fromind, int toind, int step, T val)
{
	*npval=0;
	for (int i=fromind; i<=toind; i+=step)
		if (function[i]<val) pickedindicies[(*npval) ++]=i;
}

template<class T>
inline void PickLargerThan(int *pickedindicies, int *npval, T *function, int fromind, int toind, int step, T val)
{
	*npval=0;
	for (int i=fromind; i<=toind; i+=step)
		if (function[i]>val) pickedindicies[(*npval) ++]=i;
}

template<class T>
inline void PickEqualTo(int *pickedindicies, int *npval, T *function, int fromind, int toind, int step, T val)
{
	*npval=0;
	for (int i=fromind; i<=toind; i+=step)
		if (function[i]==val) pickedindicies[(*npval) ++]=i;
}

template<class T>
inline void PickNotEqualTo(int *pickedindicies, int *npval, T *function, int fromind, int toind, int step, T val)
{
	*npval=0;
	for (int i=fromind; i<=toind; i+=step)
		if (function[i]!=val) pickedindicies[(*npval) ++]=i;
}


template<class T>
inline T Min(const T &a, const T &b) { return (a<b ? a : b); }

template<class T>
inline T Max(const T &a, const T &b) { return (a>b ? a : b); }

template<class T>
inline int Sign(const T a)
{
	if (fabs(a)<SMALL) return 0;
	if (a>0) return 1;
	return -1;
}

inline double SIGN(double a, double b) { return ((b>=0.0) ? fabs(a) : -fabs(a)); }	//a numrec function

//inline double _MAX(double a, double b) { return ((a>b) ? a : b ); }	//a numrec function

inline int GetBit(int x, int b) { return ((x >> b) & ~(~0 << 1)); }

inline void PutBit(int *x, int b) { (*x)=((*x) | ( 1 << b)); }

inline double Sqr(double t) { return t*t; }

inline void PutBitOff(int *x, int b) { (*x)=((*x) & ( ~(1 << b))); }

inline int FindLowestBit(int x)
{
	int i=0;
	if (x==0) return -1;
	while (!(GetBit(x, i))) i++;
	return i;
}


inline void Statisticsf(double *hist, int nbins, double *ave, double *spr, int no_zero)
{
	int i;
	double s, s2, dnbins;
	s=s2=0;
	dnbins=0.0;
	for (i=0; i<nbins; i++)
	{
		if (!((no_zero) && (fabs(hist[i])<SMALL)))
		{
			s+=hist[i]; s2+=hist[i]*hist[i];
			dnbins++;
		}
	}
	if (dnbins>0)
	{
		*ave=s/dnbins;
		*spr=sqrt(s2/dnbins-(*ave)*(*ave));
	}
}



class randomgenerator_class {
public:

	randomgenerator_class(bool use_time=true) {
		if (!use_time) _ranpoint=-1;
		else _ranpoint=-time(0);
		_Initialize();
	}

	int irand(int r) {
		int i=((int) (r*drand()));
		if (i>=r) i--;
		return i;
	}

	double drand();
	void SetSeed(long val) { _ranpoint=val;_Initialize(); }
	long GetSeed() const { return _ranpoint; }

private:
	long _ranpoint;
	static const int IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2;
	static const int NTAB,NDIV;
	static const double RANEPS,RNMX,AM;
	vector<long> iv;
	long iy,idum2;
	void _Initialize();

};

// NUMERICAL RECIPES ROUTINES !

void zbrak(double (*fx)(double), double x1, double x2, int n, double xb1[],
		   double xb2[], int *nb,bool output=false);
double zbrent(double (*func)(double), double x1, double x2, double tol);


// ALL THE FOLLOWING FUNCTIONS ARE FOR MANIPULATING C POINTERS
// AND INCLUDE ALLOCATION,REALLOCATION,AND SETTINGS

static const int ALLO_DEBUG=0;

template<class T>
inline void VecEq(T* v1 , const T* v2, int n) { for (int i=0; i<n; i++) v1[i]=v2[i]; }

template<class T>
inline void VecZero(T* v1, int n)
{
	if (v1==NULL) return;
	for (int i=0; i<n; i++) v1[i]=0;
}

template<class T>
inline void MatEq(T **ma1, T **ma2, const int dim1, const int dim2)
{
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++) ma1[i] [j]=ma2[i] [j];
}

template<class T>
inline void MatZero(T **ma1, const int dim1, const int dim2)
{
	if (ma1==NULL) return;
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++) ma1[i] [j]=0;
}

template<class T>
inline void LatEq(T ***la1, T ***la2, const int dim1, const int dim2, const int dim3)
{
	int i, j, k;
	if (la2==NULL)
	{
		for (i=0; i<dim1; i++)
		{
			for (j=0; j<dim2; j++) { for (k=0; k<dim3; k++) la1[i] [j] [k]=0; }
		}
	}
	else
	{
		for (i=0; i<dim1; i++)
		{
			for (j=0; j<dim2; j++) { for (k=0; k<dim3; k++) la1[i] [j] [k]=la2[i] [j] [k]; }
		}
	}
}


template<class T>
void AlloVec(T* &ve, const int dim1)
{
	if (ALLO_DEBUG>4) printf("trying to AlloVec(%d)\n", dim1);
	ve=(T *) calloc(dim1, sizeof(T));
	if (!ve) { cerr<<"allocation error in Allovec"; exit(1); }
}

template<class T>
void AlloMat(T** &m, const int dim1, const int dim2)
{
	int i;
	if (ALLO_DEBUG>3) printf("trying to AlloMat(%d,%d)\n", dim1, dim2);
	m = (T **) calloc(dim1, sizeof(T *));
	if (!m) { cerr<<" allocation error in allomat 1"; exit(1); }
	AlloVec(m[0], dim1*dim2);
	for (i=1; i<dim1; i++) { m[i] = m[i-1]+dim2; }
}

template<class T>
void AlloLat(T*** &l, const int dim1, const int dim2, const int dim3)
{
	int i;
	l=(T ***) calloc(dim1, sizeof(T **));
	if (!l) { cerr<<"allocation error in alloLat"; exit(1); }
	for (i=0; i<dim1; i++) AlloMat(*(l+i), dim2, dim3);
}

template<class T>
inline void FreeVec(T* &pi) { if (pi!=NULL) free(pi); }

template<class T>
inline void FreeMat(T** &ppi)
{
	if (ppi==NULL) return;
	free(ppi[0]);
	free(ppi);
}


#endif
