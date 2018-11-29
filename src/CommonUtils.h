#ifndef CommonUtils_h
#define CommonUtils_h

#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include "vg.pb.h"

namespace CommonUtils
{
	struct InvalidGraphException : std::runtime_error
	{
		InvalidGraphException(const char* c);
	};
	namespace inner
	{
		bool alignmentLengthCompare(const vg::Alignment* const left, const vg::Alignment* const right);
		bool alignmentScoreCompare(const vg::Alignment* const left, const vg::Alignment* const right);
		bool alignmentIncompatible(const vg::Alignment* const left, const vg::Alignment* const right);
	}
	vg::Graph LoadVGGraph(std::string filename);
	std::string ReverseComplement(std::string original);
	vg::Alignment LoadVGAlignment(std::string filename);
	std::vector<vg::Alignment> LoadVGAlignments(std::string filename);
	template <typename T, typename F>
	std::vector<T> SelectAlignments(std::vector<T> alignments, size_t maxnum, F alnGetter)
	{
		std::function<const vg::Alignment*(const T&)> f = [alnGetter](const T& aln) { return (const vg::Alignment*)alnGetter(aln); };
		std::sort(alignments.begin(), alignments.end(), [f](const T& left, const T& right) { return inner::alignmentScoreCompare(f(left), f(right)); });
		std::stable_sort(alignments.begin(), alignments.end(), [f](const T& left, const T& right) { return inner::alignmentLengthCompare(f(left), f(right)); });
		std::vector<T> result;
		assert(f(alignments[0])->sequence().size() > f(alignments.back())->sequence().size() || (f(alignments[0])->sequence().size() == f(alignments.back())->sequence().size() && f(alignments[0])->score() <= f(alignments.back())->score()));
		for (size_t i = 0; i < alignments.size(); i++)
		{
			const vg::Alignment* const aln = f(alignments[i]);
			if (!std::any_of(result.begin(), result.end(), [aln, f](const T& existing) { return inner::alignmentIncompatible(f(existing), aln); }))
			{
				result.push_back(alignments[i]);
			}
			if (result.size() == maxnum) break;
		}
		return result;
	}
	std::vector<vg::Alignment> SelectAlignments(std::vector<vg::Alignment> alns, size_t maxnum);
	std::vector<vg::Alignment*> SelectAlignments(std::vector<vg::Alignment*> alns, size_t maxnum);
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
