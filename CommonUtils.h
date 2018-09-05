#ifndef CommonUtils_h
#define CommonUtils_h

#include <string>
#include <vector>
#include <sstream>

namespace CommonUtils
{
	std::string ReverseComplement(std::string original);
}

class BufferedWriter : std::ostream
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
	static FlushClass Flush;
private:
	std::ostream* stream;
	std::stringstream stringstream;
};

#endif
