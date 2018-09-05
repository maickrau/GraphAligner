#include "CommonUtils.h"
#include "stream.hpp"

namespace CommonUtils
{

	std::string ReverseComplement(std::string str)
	{
		std::string result;
		result.reserve(str.size());
		for (int i = str.size()-1; i >= 0; i--)
		{
			switch (str[i])
			{
				case 'A':
				case 'a':
				result += 'T';
				break;
				case 'C':
				case 'c':
				result += 'G';
				break;
				case 'T':
				case 't':
				result += 'A';
				break;
				case 'G':
				case 'g':
				result += 'C';
				break;
				case 'N':
				case 'n':
				result += 'N';
				break;
				case 'U':
				case 'u':
				result += 'A';
				break;
				case 'R':
				case 'r':
				result += 'Y';
				break;
				case 'Y':
				case 'y':
				result += 'R';
				break;
				case 'K':
				case 'k':
				result += 'M';
				break;
				case 'M':
				case 'm':
				result += 'K';
				break;
				case 'S':
				case 's':
				result += 'S';
				break;
				case 'W':
				case 'w':
				result += 'W';
				break;
				case 'B':
				case 'b':
				result += 'V';
				break;
				case 'V':
				case 'v':
				result += 'B';
				break;
				case 'D':
				case 'd':
				result += 'H';
				break;
				case 'H':
				case 'h':
				result += 'D';
				default:
				assert(false);
			}
		}
		return result;
	}

}

BufferedWriter::BufferedWriter() : stream(nullptr) {};
BufferedWriter::BufferedWriter(std::ostream& stream) : stream(&stream) {};
BufferedWriter& BufferedWriter::operator<<(FlushClass f)
{
	if (stream == nullptr) return *this;
	flush();
	return *this;
}
void BufferedWriter::flush()
{
	if (stream == nullptr) return;
	stringstream << std::endl;
	(*stream) << stringstream.str();
	stringstream.str("");
}
