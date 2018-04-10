#ifndef CommonUtils_h
#define CommonUtils_h

#include <string>
#include <vector>
#include <sstream>
#include "vg.pb.h"

namespace CommonUtils
{
	vg::Graph LoadVGGraph(std::string filename);
	std::string ReverseComplement(std::string original);
	vg::Alignment LoadVGAlignment(std::string filename);
	std::vector<vg::Alignment> LoadVGAlignments(std::string filename);
}

class BufferedWriter : std::ostream
{
public:
	class FlushClass {};
	BufferedWriter();
	BufferedWriter(std::ostream& stream);
	template <typename T>
	BufferedWriter& operator<<(T obj)
	{
		if (stream == nullptr) return *this;
		stringstream << obj;
		return *this;
	}
	BufferedWriter& operator<<(FlushClass f);
	void flush();
	static FlushClass Flush;
private:
	std::ostream* stream;
	std::stringstream stringstream;
};

#endif
