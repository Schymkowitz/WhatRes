// SequenceAnalysis, Copyright Switch Lab 2020
// Joost Schymkowitz

#ifndef NUMHIST_H
#define NUMHIST_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <utility>
#include <math.h>
#include <algorithm>

using namespace std;

class NumberHistogram
{
	vector<double> counts;
	vector<double> fractions;
	int bins;
	double binsize;
	double min_val, max_val;
	double totalcount;
	
public:
	NumberHistogram(int number_of_bins, double _min, double _max)
		:bins(number_of_bins), min_val(_min), max_val(_max), totalcount(0.), binsize((_max-_min)/number_of_bins)
	     ,counts(number_of_bins, 0.), fractions(number_of_bins, 0.)
	{
//		cout << "histogram created, bin number = " << bins << " " << counts.size() << endl;
	}
	
	int findbin(double num)
	{
		if(num < min_val)
		{
			cout << "number below range in histogram " << num << endl;
			return 0;
		}
		else if(num > max_val)
		{
			cout << "number over range in histogram " << num << endl;
			return bins;
		}
		int temp = floor((num - min_val) / binsize);
//		cout << "bin: " << temp << endl;
		return temp;
	}
	void Add(double num) // push a number into the histogram
	{
//		cout << "push into histogram" << endl;
		int bin = findbin(num);
//		cout << "accessing counts " << endl;
		counts[bin]++;
		totalcount++;
	}
	void Totalise()  // convert counts to fractions
	{
		for(int i = 0; i < bins; i++)
		{
			fractions[i] = counts[i] / totalcount;
		}
	}
	
	void Normalise() // specifically relevant to volume based binning
	{
		for(int i= 0; i < bins; i++)
		{
			double midbin = i * binsize + 0.5 * binsize;
			fractions[i] = fractions[i] / pow(midbin, 3);
		}
	}
	
	void PrintFractions(ofstream& outer)
	{
		for(int i=0; i< bins; i++)
		{
			outer << fractions[i] << "\t";
		}
	}
	void PrintCounts(ofstream& outer)
	{
		for(int i=0; i < bins; i++)
		{
			outer << counts[i] << "\t";
		}
		outer << totalcount;
	}
	double Get_totalcount() {return totalcount;}
};

#endif

