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
	namespace inner
	{
		//an overlap which is larger than the fraction cutoff of the smaller alignment means the alignments are incompatible
		//eg alignments 12000bp and 15000bp, overlap of 12000*0.05 = 600bp means they are incompatible
		constexpr float OverlapIncompatibleFractionCutoff = 0.05;
		template <typename T>
		bool alignmentCompare(const T& left, const T& right)
		{
			if (left.alignmentEnd - left.alignmentStart < right.alignmentEnd - right.alignmentStart) return true;
			if (right.alignmentEnd - right.alignmentStart < left.alignmentEnd - left.alignmentStart) return false;
			if (left.alignmentScore < right.alignmentScore) return true;
			return false;
		}
		template <typename T>
		bool alignmentIncompatible(const T& left, const T& right)
		{
			auto minOverlapLen = std::min(left.alignmentEnd - left.alignmentStart, right.alignmentEnd - right.alignmentStart) * OverlapIncompatibleFractionCutoff;
			size_t leftStart = left.alignmentStart;
			size_t leftEnd = left.alignmentEnd;
			size_t rightStart = right.alignmentStart;
			size_t rightEnd = right.alignmentEnd;
			if (leftStart > rightStart)
			{
				std::swap(leftStart, rightStart);
				std::swap(leftEnd, rightEnd);
			}
			int overlap = 0;
			assert(leftStart <= rightStart);
			if (leftEnd > rightStart) overlap = leftEnd - rightStart;
			return overlap > minOverlapLen;
		}
	}
	vg::Graph LoadVGGraph(std::string filename);
	char Complement(char original);
	std::string ReverseComplement(std::string original);
	vg::Alignment LoadVGAlignment(std::string filename);
	std::vector<vg::Alignment> LoadVGAlignments(std::string filename);
	template <typename T>
	std::vector<T> SelectAlignments(std::vector<T> alignments, size_t maxnum)
	{
		std::sort(alignments.begin(), alignments.end(), [](const T& left, const T& right) { return inner::alignmentCompare(left, right); });
		std::vector<T> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			if (!std::any_of(result.begin(), result.end(), [&alignments, i](const T& existing) { return inner::alignmentIncompatible(existing, alignments[i]); }))
			{
				result.emplace_back(std::move(alignments[i]));
			}
			if (result.size() == maxnum) break;
		}
		return result;
	}
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
	bool inputDiscarded() const;
	static FlushClass Flush;
private:
	std::ostream* stream;
	std::stringstream stringstream;
};

#endif
