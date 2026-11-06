// commands.h
// SequenceAnalysis, Copyright Switch Lab 2010
// Joost Schymkowitz
//
// Base class for all analysis commands.
//
// ---------------------------------------------------------------------
// Architecture overview
// ---------------------------------------------------------------------
// The framework uses a Command Pattern:
//
//   command_switch             holds a registry of all commands
//   sequence_command (this)    abstract base every command inherits from
//   <specific>_command         concrete implementation of one analysis
//
// Each concrete command:
//   1. Inherits from sequence_command
//   2. Passes a label string and explanation to the base constructor
//   3. Implements Execute() with the actual analysis logic
//   4. Optionally overrides SetOptionsRequired() to document which
//      Options.txt fields it reads
//
// Commands are invoked from the command line as:
//   whatres-features -do <label> -in <file> -opt <optionsfile>
//
// The label string passed to the base constructor must match exactly
// what the user types after -do.
//
// ---------------------------------------------------------------------
// Writing a new command: minimal example
// ---------------------------------------------------------------------
//
//   // my_command.h
//   #include "commands.h"
//   class my_command : public sequence_command {
//   public:
//       my_command(Option_class& O)
//           : sequence_command("my_command", "one-line description", O) {}
//       void Execute() override;
//   };
//
//   // my_command.cpp
//   #include "my_command.h"
//   void my_command::Execute() { /* analysis logic here */ }
//
// Then in command_switch.cpp:
//   #include "my_command.h"
//   _commands.emplace_back(CMD(my_command));   // in the constructor body
// ---------------------------------------------------------------------

#ifndef COMMANDS_H
#define COMMANDS_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

#include "optionclass.h"

using namespace std;

// sequence_command: abstract base class for all analysis commands.
//
// Provides:
//   - A label (used for command-line dispatch and GUI display)
//   - An explanation string (shown in help / GUI)
//   - An Option_class instance (populated from Options.txt before Execute())
//   - Standard I/O helpers (PrintInputFiles, ReplaceInputFiles, etc.)
//   - A list of option names the command requires (for documentation)
class sequence_command
{

protected:
	string label;               // Command name used for -do dispatch (e.g. "tango")
	string explanation;         // One-line description shown in GUI and -help output
	Option_class Option;        // Configuration populated by Read_Options() before Execute()
	vector<string> options_required; // Names of Options.txt fields this command reads
	bool include_in_gui;        // If true, this command is exposed in the web GUI

public:
	// Construct with a fixed label, explanation, and initial option set.
	// Calls SetOptionsRequired() so subclasses can populate options_required
	// at construction time via an override.
	// By default commands are NOT included in the GUI; subclasses opt in by
	// setting include_in_gui = true in their constructor body.
	sequence_command(const string& l, const string& e, Option_class &O)
		: label(l), explanation(e), Option(O), include_in_gui(false)
	{
		SetOptionsRequired();
	};

	// Virtual destructor, required because command_switch deletes commands
	// through a base-class pointer.
	virtual ~sequence_command() {}

	// ---------------------------------------------------------------------------
	// Interface every command must implement
	// ---------------------------------------------------------------------------

	// Run the analysis. Called after Read_Options() and ReplaceInputFiles().
	// Input files are available in Option.inputfiles.
	virtual void Execute() = 0;

	// ---------------------------------------------------------------------------
	// Optional overrides
	// ---------------------------------------------------------------------------

	// Override to return a detailed, multi-paragraph help string for this
	// command: what it does, the input files it expects, what it writes,
	// the scientific method / reference.  Both the `-help` command-line
	// output and the generated DOCUMENTATION.md are built from this, so it
	// is the single source of truth for a command's documentation.
	// Empty by default (the command then shows only its one-line label).
	virtual string HelpText() { return ""; }

	// Print formatted help for this command to stdout: the dispatch label,
	// the one-line explanation, the detailed HelpText(), and the
	// Options.txt fields the command reads.  Invoked by `-do <cmd> -help`.
	virtual void GetHelp()
	{
		cout << "\n=== " << label << " ===\n"
		     << explanation << "\n";
		string ht = HelpText();
		if(!ht.empty()) cout << "\n" << ht << "\n";
		if(!options_required.empty())
		{
			cout << "\nOptions.txt fields read by this command:\n";
			for(size_t i = 0; i < options_required.size(); ++i)
				cout << "  " << options_required[i] << "\n";
		}
		cout << endl;
	}

	// Return false if required input files are missing or malformed.
	// Called before Execute() to give a clean error message.
	virtual bool checkInputFiles() { return true; }

	// Override to give a longer description (currently unused by dispatcher).
	virtual void SetExplanation() {}

	// Override to push the names of Options.txt fields this command reads
	// into options_required. Called from the constructor.
	// Example:
	//   void SetOptionsRequired() override {
	//       options_required.push_back("tango_window_tresh");
	//   }
	virtual void SetOptionsRequired() {}

	// ---------------------------------------------------------------------------
	// Provided utilities (do not override)
	// ---------------------------------------------------------------------------

	// Execute once per input file listed in Option.inputfiles.
	void Run()
	{
		int sizer = Option.inputfiles.size();
		for(int i = 0; i < sizer; i++)
		{
			Execute();
		}
	}

	// Return the dispatch label (the string used after -do on the command line).
	string GetLabel() { return label; }
	// Match the -do dispatch name. Override to accept aliases (e.g. an old name).
	virtual bool MatchesLabel(const string& c) { return label == c; }

	// Return the one-line explanation string.
	string GetExplanation() { return explanation; }

	// Return whether this command should appear in the web GUI.
	bool GetIncludeInGUI() { return include_in_gui; }

	// Write the current option values to PSX_usedOptions.txt for reproducibility.
	void PrintOptions()
	{
		ofstream o("PSX_usedOptions.txt");
		Option.PrintOptions(o);
		o.close();
	}

	// Populate Option from an XML options file.
	// Returns false if the file cannot be opened (uses built-in defaults then).
	bool Read_Options(const string& FN)
	{
		return Option.Read_Options(FN);
	}

	// Print the list of input files currently registered in Option to stdout.
	void PrintInputFiles()
	{
		Option.PrintInputFiles();
	}

	// Replace the input file list in Option with the supplied vector.
	// Called by command_switch when -in arguments are given on the command line.
	void ReplaceInputFiles(vector<string>& newfiles)
	{
		Option.ReplaceInputFiles(newfiles);
	}

	// Replace the mutation list in Option.
	// Called by command_switch when -mut arguments are given on the command line.
	void ReplaceMutations(vector<string>& muts)
	{
		Option.ReplaceMutations(muts);
	}

	// Replace the output prefix string (default "PSX") in Option.
	// Called by command_switch when -pre is given on the command line.
	void ReplacePrefix(const string& pre)
	{
		Option.ReplacePrefix(pre);
	}

	// Return a stringstream listing all option names this command requires,
	// formatted as <option>name\n lines for GUI / help output.
	stringstream GetOptionsRequired()
	{
		int optsize = options_required.size();
		stringstream ss;
		for(int i = 0; i < optsize; i++)
		{
			ss << "<option>" << options_required[i] << "\n";
		}
		return ss;
	}
};

#endif
