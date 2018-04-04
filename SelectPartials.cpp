#include <algorithm>
#include <string>
#include <unordered_map>
#include "CommonUtils.h"
#include "stream.hpp"

//an overlap which is larger than the fraction cutoff of the smaller alignment means the alignments are incompatible
//eg alignments 12000bp and 15000bp, overlap of 12000*0.05 = 600bp means they are incompatible
const float OverlapIncompatibleFractionCutoff = 0.05;

std::unordered_map<std::string, std::vector<vg::Alignment>> splitAlignments(const std::vector<vg::Alignment>& alignments)
{
	std::unordered_map<std::string, std::vector<vg::Alignment>> result;
	for (auto alignment : alignments)
	{
		result[alignment.name()].push_back(alignment);
	}
	return result;
}

//longer alignments are better
bool alignmentLengthCompare(const vg::Alignment& left, const vg::Alignment& right)
{
	return left.sequence().size() > right.sequence().size();
}

//lower scores are better
bool alignmentScoreCompare(const vg::Alignment& left, const vg::Alignment& right)
{
	return left.score() < right.score();
}

bool alignmentIncompatible(const vg::Alignment& left, const vg::Alignment& right)
{
	auto minOverlapLen = std::min(left.sequence().size(), right.sequence().size()) * OverlapIncompatibleFractionCutoff;
	auto leftStart = left.query_position();
	auto leftEnd = leftStart + left.sequence().size();
	auto rightStart = right.query_position();
	auto rightEnd = rightStart + right.sequence().size();
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

std::vector<vg::Alignment> selectAlignments(std::vector<vg::Alignment> alignments)
{
	std::sort(alignments.begin(), alignments.end(), alignmentScoreCompare);
	std::stable_sort(alignments.begin(), alignments.end(), alignmentLengthCompare);
	std::vector<vg::Alignment> result;
	assert(alignments[0].sequence().size() > alignments.back().sequence().size() || (alignments[0].sequence().size() == alignments.back().sequence().size() && alignments[0].score() <= alignments.back().score()));
	for (size_t i = 0; i < alignments.size(); i++)
	{
		vg::Alignment& aln = alignments[i];
		if (!std::any_of(result.begin(), result.end(), [&aln](const vg::Alignment& existing) { return alignmentIncompatible(existing, aln); }))
		{
			result.push_back(aln);
		}
	}
	return result;
}

int main(int argc, char** argv)
{
	std::string infilename {argv[1]};
	std::string outfilename {argv[2]};

	auto alignments = CommonUtils::LoadVGAlignments(infilename);
	std::unordered_map<std::string, std::vector<vg::Alignment>> alignmentsPerRead = splitAlignments(alignments);
	std::vector<vg::Alignment> totalSelected;
	for (auto pair : alignmentsPerRead)
	{
		std::vector<vg::Alignment> selected = selectAlignments(pair.second);
		totalSelected.insert(totalSelected.end(), selected.begin(), selected.end());
	}
	std::ofstream outfile { outfilename,  std::ios::out | std::ios::binary };
	stream::write_buffered(outfile, totalSelected, 0);
}
