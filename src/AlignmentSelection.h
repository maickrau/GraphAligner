#ifndef AlignmentSelection_h
#define AlignmentSelection_h

#include <cmath>
#include <cstdint>
#include <vector>
#include <functional>
#include "vg.pb.h"
#include "GraphAlignerCommon.h"
#include "EValue.h"

namespace AlignmentSelection
{
	extern double OverlapIncompatibleFractionCutoff;
	struct SelectionOptions
	{
		size_t graphSize;
		size_t readSize;
		double ECutoff;
		EValueCalculator EValueCalc;
		double AlignmentScoreFractionCutoff;
		int minAlignmentScore;
	};
	void RemoveDuplicateAlignments(const AlignmentGraph& graph, std::vector<AlignmentResult::AlignmentItem>& alignments);
	std::vector<AlignmentResult::AlignmentItem> SelectAlignments(const std::vector<AlignmentResult::AlignmentItem>& alignments, SelectionOptions options);
	bool alignmentIncompatible(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right);

	std::vector<AlignmentResult::AlignmentItem> SelectAlignmentScore(const std::vector<AlignmentResult::AlignmentItem>& alignments, double score, const EValueCalculator& EValueCalc);
	std::vector<AlignmentResult::AlignmentItem> SelectECutoff(const std::vector<AlignmentResult::AlignmentItem>& alignments, size_t m, size_t n, double cutoff, const EValueCalculator& EValueCalc);
	std::vector<AlignmentResult::AlignmentItem> SelectAlignmentFractionCutoff(const std::vector<AlignmentResult::AlignmentItem>& alignments, double cutoff, const EValueCalculator& EValueCalc);
	void AddMappingQualities(std::vector<AlignmentResult::AlignmentItem>& alignments);

};

#endif
