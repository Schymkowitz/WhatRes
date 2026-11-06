// command_switch.cpp
// SequenceAnalysis, Copyright Switch Lab 2010
// Joost Schymkowitz

#include "command_switch.h"

#include "optionclass.h"
#include "pdb_quality_check_command.h"
#include "pdb_surface_analysis_command.h"
#include "PDB_dist_command.h"
#include "structure_survey_command.h"

// Convenience macro: wraps a raw new in a unique_ptr<sequence_command>.
#define CMD(T) unique_ptr<sequence_command>(new T(O))

command_switch::~command_switch() = default;

command_switch::command_switch()
{
	Option_class O;

	_commands.emplace_back(CMD(pdb_quality_check_command));
	_commands.emplace_back(CMD(pdb_surface_analysis_command));
	_commands.emplace_back(CMD(PDB_dist_command));
	_commands.emplace_back(CMD(structure_survey_command));
}

#undef CMD

void command_switch::RunCommandLine(int my_argc, char* my_argv[])
{
	int sizer = _commands.size(), choice=-1;
	vector<string> helplist, inputs, muts;
	bool command=false, inputf=false, optf=false, mutations=false, pref=false;
	string opt = "Options.txt", com, pre;

	for(int i=0; i<my_argc; i++) helplist.push_back(string(my_argv[i]));

	// `-help` / `--help` anywhere on the command line requests help rather
	// than execution: `-do <cmd> -help` prints that command's help, a bare
	// `-help` lists every command.
	bool helpRequested = false;
	for(size_t i=0; i<helplist.size(); ++i)
		if(helplist[i]=="-help" || helplist[i]=="--help") helpRequested=true;

	for(int i=0; i<my_argc; i++)
	{
		if(helplist[i] != "-do" && helplist[i] != "-in" && helplist[i] != "-opt" && helplist[i] != "-mut" && helplist[i] != "-pre") continue;
		i++;
		if(i==my_argc)
		{
			cout << "missing arguments, aborting command line operation" << endl;
			return;
		}

		if(helplist[i-1] == "-do")
		{
			command=true;
			com = helplist[i];
		}
		else if(helplist[i-1] == "-in")
		{
			inputf=true;
			inputs.push_back(helplist[i]);
			int e = i+1;
			while(e<my_argc && helplist[e][0]!='-')
			{
				inputs.push_back(helplist[e]);
				e++;
			}
		}
		else if(helplist[i-1] == "-opt")
		{
			optf=true;
			opt = helplist[i];
		}
		else if(helplist[i-1] == "-pre")
		{
			pref = true;
			pre = helplist[i];
		}
		else if(helplist[i-1] == "-mut")
		{
			mutations = true;
			muts.push_back(helplist[i]);
			int e = i+1;
			while(e<my_argc && helplist[e][0]!='-')
			{
				muts.push_back(helplist[e]);
				e++;
			}
		}
	}
	cout<<"listsize= "<<muts.size() << endl;
	for(int i=0; i<sizer; i++)
	{
		if(_commands[i]->MatchesLabel(com))
		{
			choice = i;
			cout<<"label: "<<com<<" found."<<endl;
			break;
		}
	}

	// Help request: print and stop, before any options / input handling.
	if(helpRequested)
	{
		if(choice>=0) _commands[choice]->GetHelp();
		else          PrintCommandList();
		return;
	}

	if(choice==-1)
	{
		cout<<"command line argument not understood "<<com << endl;
		return;
	}

	if(!_commands[choice]->Read_Options(opt))
	{
		cout << "Option file not found: " << opt << " , aborting." << endl;
		return;
	}

	if(!inputf)
	{
		cout << "no inputfiles specified, reading input from option file."<<endl;
	}
	else
	{
		_commands[choice]->ReplaceInputFiles(inputs);
	}

	if(mutations)
	{
		_commands[choice]->ReplaceMutations(muts);
	}

	if(pref)
	{
		_commands[choice]->ReplacePrefix(pre);
	}

	cout<<"using optionsfile: "<<opt<<endl;
	_commands[choice]->PrintInputFiles();
	_commands[choice]->PrintOptions();
	_commands[choice]->Execute();

}

void command_switch::PrintCommandList()
{
	cout << "\nwhatres-features: " << _commands.size() << " commands\n";
	cout << "Run  whatres-features -do <command> -help  for detail on one.\n\n";
	for(size_t i=0; i<_commands.size(); ++i)
		cout << "  " << _commands[i]->GetLabel() << "\n"
		     << "      " << _commands[i]->GetExplanation() << "\n";
	cout << endl;
}
