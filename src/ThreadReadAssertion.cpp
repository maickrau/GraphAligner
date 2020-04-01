#include <iostream>
#include <sstream>
#include <string_view>
#include "ThreadReadAssertion.h"

namespace ThreadReadAssertion
{
	thread_local int currentnodeID;
	thread_local bool currentreverse;
	thread_local size_t currentseqPos;
	thread_local size_t currentmatchLen;
	thread_local size_t currentnodeOffset;
	thread_local std::string_view currentRead;
	void signal(int signal)
	{
		std::stringstream msg;
		msg << "Signal " << signal << ". Read: " << currentRead << ". Seed: " << assertGetSeedInfo();
		std::cerr << msg.str() << std::endl;
		std::abort();
	}
	void setRead(const std::string& readName)
	{
		currentRead = std::string_view(readName.data(), readName.size());
	}
	void setSeed(int nodeID, bool reverse, size_t seqPos, size_t matchLen, size_t nodeOffset)
	{
		currentnodeID = nodeID;
		currentreverse = reverse;
		currentseqPos = seqPos;
		currentmatchLen = matchLen;
		currentnodeOffset = nodeOffset;
	}
	void assertFailed(const char* expression, const char* file, int line)
	{
		std::stringstream msg;
		msg << file << ":" << line << ": Assertion '" << expression << "' failed. Read: " << currentRead << ". Seed: " << assertGetSeedInfo();
		std::cerr << msg.str() << std::endl;
		throw AssertionFailure {};
	}	
	std::string assertGetSeedInfo()
	{
		return std::to_string(currentnodeID) + (currentreverse ? "-" : "+") + "," + std::to_string(currentseqPos) + "," + std::to_string(currentmatchLen) + "," + std::to_string(currentnodeOffset);
	}
}
