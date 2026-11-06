// main.cpp
// SequenceAnalysis, Copyright Switch Lab 2010
// Joost Schymkowitz
//
//   whatres-features -do <command> -in <pdb_list.txt> -opt <Options.txt>
//   whatres-features -help

#include <iostream>
#include <stdexcept>
#include <string>

#include "command_switch.h"

using namespace std;

int main(int argc, char* argv[])
{
	try
	{
		command_switch CS;
		if(argc >= 3 && string(argv[1]) == "-do")
		{
			CS.RunCommandLine(argc, argv);
			return 0;
		}
		CS.PrintCommandList();
		if(argc == 1) return 0;
		const string a(argv[1]);
		return (a == "-help" || a == "--help" || a == "help") ? 0 : 1;
	}
	catch(const exception& e)
	{
		cerr << "Error: " << e.what() << endl;
		return 1;
	}
}
