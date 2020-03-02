#ifndef ThreadReadAssertion_h
#define ThreadReadAssertion_h

#include <string>

namespace ThreadReadAssertion
{
	class AssertionFailure
	{
	};
	void setRead(const std::string& readName);
	void setSeed(int nodeID, bool reverse, size_t seqPos, size_t matchLen, size_t nodeOffset);
	void assertFailed(const char* expression, const char* file, int line);
	void signal(int signal);
	std::string assertGetSeedInfo();
}

#endif

#ifdef assert
#undef assert
#endif

#ifndef NDEBUG

//https://stackoverflow.com/questions/9701229/c-assert-implementation-in-assert-h
#define assert(expression) (void)((expression) || (ThreadReadAssertion::assertFailed(#expression, __FILE__, __LINE__),0))
#define assertSetRead(name, nodeid, reverse, seqpos, matchlen, nodeoffset) { ThreadReadAssertion::setRead(name); ThreadReadAssertion::setSeed(nodeid, reverse, seqpos, matchlen, nodeoffset); }
#define assertSetNoRead(name) { ThreadReadAssertion::setRead(name); ThreadReadAssertion::setSeed(0, 0, 0, 0, 0); }

#else

#define assert(ignore) ((void)0)
#define assertSetRead(name, nodeid, reverse, seqpos, matchlen, nodeoffset) ((void)0)
#define assertSetNoRead(name) ((void)0)

#endif
