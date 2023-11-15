#ifndef CommonUtils_h
#define CommonUtils_h

#include <algorithm>
#include <functional>
#include <string>
#include <vector>
#include <sstream>
#include "vg.pb.h"

namespace CommonUtils
{
	struct InvalidGraphException : std::runtime_error
	{
		InvalidGraphException(const char* c);
		InvalidGraphException(std::string c);
	};
	vg::Graph LoadVGGraph(std::string filename);
	char Complement(char original);
	std::string ReverseComplement(std::string original);
	vg::Alignment LoadVGAlignment(std::string filename);
	std::vector<vg::Alignment> LoadVGAlignments(std::string filename);
}

class BufferedWriter
{
public:
	class FlushClass {};
	BufferedWriter();
	BufferedWriter(std::ostream& stream);
	BufferedWriter(const BufferedWriter& other) = default;
	BufferedWriter(BufferedWriter&& other) = default;
	BufferedWriter& operator=(const BufferedWriter& other) = default;
	BufferedWriter& operator=(BufferedWriter&& other) = default;
	template <typename T>
	BufferedWriter& operator<<(T obj)
	{
		if (stream == nullptr) return *this;
		stringstream << obj;
		return *this;
	}
	BufferedWriter& operator<<(FlushClass f);
	void flush();
	bool inputDiscarded() const;
	static FlushClass Flush;
private:
	std::ostream* stream;
	std::stringstream stringstream;
};

#if !defined(__clang__)
namespace std
{
	template <> struct hash<__uint128_t>
	{
		size_t operator()(__uint128_t x) const
		{
			return (size_t)x ^ (size_t)(x >> 64);
		}
	};
}
#endif

#endif
