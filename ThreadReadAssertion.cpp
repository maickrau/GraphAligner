#include <iostream>
#include <sstream>
#include "ThreadReadAssertion.h"

namespace ThreadReadAssertion
{
	thread_local std::string currentRead;
	void setRead(const std::string& readName)
	{
		currentRead = readName;
	}
	void assertFailed(const char* expression, const char* file, int line)
	{
		std::stringstream msg;
		msg << file << ":" << line << ": Assertion '" << expression << "' failed. Read: " << currentRead;
		std::cerr << msg.str() << std::endl;
		throw AssertionFailure {};
	}	
}
