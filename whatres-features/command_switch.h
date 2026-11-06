// command_switch.h
// SequenceAnalysis, Copyright Switch Lab 2010
// Joost Schymkowitz
//
// Command dispatcher: holds the registered commands and runs one of them
// from the command line:
//
//   whatres-features -do <command> -in <file1> [file2 ...] [-opt options.txt] [-pre <prefix>]
//
//   -do   name of the command to run (must match a registered label)
//   -in   input file(s); for these commands a text file listing one PDB path per line
//   -opt  path to the options file (default: Options.txt)
//   -pre  output file prefix (default: PSX)
//   -help with -do: help for that command; alone: the list of commands
#ifndef COMMAND_SWITCH_H
#define COMMAND_SWITCH_H

#include <iostream>
#include <string>
#include <memory>
#include <vector>

#include "commands.h"

using namespace std;

class command_switch
{
public:
	command_switch();    // registers every command
	~command_switch();

	// Parse argc/argv for -do / -in / -opt / -mut / -pre, find the command, run it.
	void RunCommandLine(int my_argc, char* my_argv[]);

	// Print every command's label and one-line explanation.
	void PrintCommandList();

private:
	vector<unique_ptr<sequence_command>> _commands;
};

#endif
