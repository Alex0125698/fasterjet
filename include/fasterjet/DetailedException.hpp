/**
 * \file DetailedException.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include <string>
#include <cassert>

class DetailedException
{
private:
	const int line;
	const std::string file;
	const std::string message;
public:
	DetailedException(std::string details, const int line, const char* file)
		: line(line), file(file), message(details.c_str())
	{
	}
	int getLine() { return line; }
	std::string getFile() { return file; }
	std::string getMessage() { return message; }
	std::string getDetails()
	{
		return "ERROR: " + message + " | " + file.substr(file.find_last_of('\\') + 1) + ":" + std::to_string(line);
	}
};

// use this to automatically add line + file details
#define DETAILEDEXCEPTION(x) DetailedException(x, __LINE__, __FILE__)
