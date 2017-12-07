#ifndef ThreadReadAssertion_h
#define ThreadReadAssertion_h

#include <string>

namespace ThreadReadAssertion
{
	class AssertionFailure
	{
	};
	void setRead(const std::string& readName);
	void assertFailed(const char* expression, const char* file, int line);
	void signal(int signal);
}

#endif

#ifndef NDEBUG

#ifdef assert
#undef assert
#endif

//https://stackoverflow.com/questions/9701229/c-assert-implementation-in-assert-h
#define assert(expression) (void)((expression) || (ThreadReadAssertion::assertFailed(#expression, __FILE__, __LINE__),0))
#define assertSetRead(name) ThreadReadAssertion::setRead(name)

#else

#define assert(ignore) ((void)0)
#define assertSetRead(ignore) ((void)0)

#endif
