// LineReader.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// GetLine: read one line of text from a stream. A line ends at the first '\r'
// or '\n', which is consumed and not stored, so a Windows "\r\n" ending yields
// an extra empty line. Used by the PDB reader, the PDB list readers and the
// options-file parser.
#pragma once

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

inline std::istream& GetLine(std::istream& in, std::string& line)
{
	line.clear();
	char c;
	while (in.get(c))
	{
		if (c == '\r' || c == '\n') break;
		line += c;
	}
	return in;
}
