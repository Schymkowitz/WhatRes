// OptionFile.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Reader for options files: one "<tag>value;" entry per line.
//   - a line that contains "//" anywhere is a comment and is skipped
//   - a line without '<', '>' and ';' is skipped
//   - the tag is the text between the first '<' and the first '>'
//   - the value is the text between '>' and the first ';'. If the line holds a
//     comma, the value is split on commas into several values.
//   - a value "#" means "keep the default"
// The first entry with a given tag wins, except in GetAll, which collects the
// values of every entry with that tag.
#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "LineReader.h"

class OptionFile
{
public:
	// Parse a file. Returns false if it cannot be opened.
	bool Read(const std::string& path)
	{
		std::ifstream in(path.c_str());
		if (!in) return false;
		std::string line;
		while (!in.eof())
		{
			GetLine(in, line);
			ParseLine(line);
		}
		return true;
	}

	// Text value, surrounding whitespace removed. Unchanged if the tag is
	// absent, the value is "#", or the value is only whitespace.
	void Get(const std::string& tag, std::string& v) const
	{
		const std::string* raw = First(tag);
		if (!raw || *raw == "#") return;
		const char* ws = " \t\r\n";
		const size_t b = raw->find_first_not_of(ws);
		if (b == std::string::npos) return;
		v = raw->substr(b, raw->find_last_not_of(ws) - b + 1);
	}

	// true / TRUE / false / FALSE; any other text leaves v unchanged.
	void Get(const std::string& tag, bool& v) const
	{
		const std::string* raw = First(tag);
		if (!raw) return;
		if (*raw == "true" || *raw == "TRUE") v = true;
		else if (*raw == "false" || *raw == "FALSE") v = false;
	}

	// Numbers, read as a stream reads them (an unreadable value gives 0).
	template <class T>
	void Get(const std::string& tag, T& v) const
	{
		const std::string* raw = First(tag);
		if (!raw || *raw == "#") return;
		std::istringstream s(*raw);
		s >> v;
	}

	// Every value of every entry with this tag, "#" values left out.
	void GetAll(const std::string& tag, std::vector<std::string>& v) const
	{
		v.clear();
		bool seen = false;
		for (const Entry& e : entries_)
		{
			if (e.tag != tag) continue;
			seen = true;
			for (const std::string& x : e.values)
				if (x != "#") v.push_back(x);
		}
		if (!seen) std::cout << "label " << tag << " not found" << std::endl;
	}

private:
	struct Entry
	{
		std::string tag;
		std::vector<std::string> values;
	};
	std::vector<Entry> entries_;

	const std::string* First(const std::string& tag) const
	{
		for (const Entry& e : entries_)
			if (e.tag == tag && !e.values.empty()) return &e.values[0];
		return nullptr;
	}

	void ParseLine(const std::string& line)
	{
		if (line.find("//") != std::string::npos) return;
		const size_t lt = line.find('<'), gt = line.find('>'), semi = line.find(';');
		if (lt == std::string::npos || gt == std::string::npos || semi == std::string::npos) return;

		Entry e;
		e.tag = line.substr(lt + 1, gt - lt - 1);
		size_t from = gt + 1;                       // start of the current value
		size_t comma = line.find(',');
		while (comma != std::string::npos)
		{
			e.values.push_back(line.substr(from, comma - from));
			from = comma + 1;
			comma = line.find(',', from);
		}
		e.values.push_back(line.substr(from, semi - from));
		entries_.push_back(e);
	}
};
