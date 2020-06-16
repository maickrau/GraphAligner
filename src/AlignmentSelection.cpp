#include <algorithm>
#include <cstdint>
#include <cmath>
#include "AlignmentSelection.h"

//an overlap which is larger than the fraction cutoff of the smaller alignment means the alignments are incompatible
//eg alignments 12000bp and 15000bp, overlap of 12000*0.05 = 600bp means they are incompatible
const float OverlapIncompatibleFractionCutoff = 0.05;

//model the alignment as one sequence with the alphabet {match, mismatch}
//with random alignments having P(match) = P(mismatch) = 0.5
//and match having score +1 and mismatch having score -2
//then use Karlin-Altschul equation to calculate E
const float lambda = 0.48121182506; //ln(1+sqrt(5)) - ln(2)
const float K = 29.8318175495*1.25982891379;

namespace AlignmentSelection
{
	//approx ~ n(matches) - 2 * n(mismatches), higher is better
	double S(const AlignmentResult::AlignmentItem& aln)
	{
		return (aln.alignmentEnd - aln.alignmentStart) - aln.alignmentScore * 3;
	}

	double Evalue(size_t m, size_t n, const AlignmentResult::AlignmentItem& aln)
	{
		return K * m * n * exp(-lambda * S(aln));
	}
	bool alignmentIncompatible(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right)
	{
		auto minOverlapLen = std::min((left.alignmentEnd - left.alignmentStart), (right.alignmentEnd - right.alignmentStart)) * OverlapIncompatibleFractionCutoff;
		assert(left.alignmentStart >= 0);
		assert(right.alignmentStart >= 0);
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

	//lower E-value is better
	bool alignmentECompare(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right, size_t m, size_t n)
	{
		return Evalue(m, n, left) < Evalue(m, n, right);
	}

	bool alignmentScoreCompare(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right)
	{
		return S(left) > S(right);
	}

	//longer is better, after that lower score is better
	bool alignmentLengthCompare(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right)
	{
		if ((left.alignmentEnd - left.alignmentStart) > (right.alignmentEnd - right.alignmentStart)) return true;
		if ((right.alignmentEnd - right.alignmentStart) > (left.alignmentEnd - left.alignmentStart)) return false;
		if (left.alignmentScore < right.alignmentScore) return true;
		return false;
	}

	std::vector<AlignmentResult::AlignmentItem> SelectAlignments(const std::vector<AlignmentResult::AlignmentItem>& allAlignments, SelectionOptions options)
	{
		// roundabout to fit the signature of const ref while allowing filtering
		std::vector<AlignmentResult::AlignmentItem> filteredByE;
		if (options.ECutoff != -1)
		{
			filteredByE = SelectECutoff(allAlignments, options.graphSize, options.readSize, options.ECutoff);
		}
		const std::vector<AlignmentResult::AlignmentItem>& alignments { (options.ECutoff != -1) ? filteredByE : allAlignments };
		switch(options.method)
		{
			case GreedyLength:
				return GreedySelectAlignments(alignments, alignmentLengthCompare);
			case GreedyScore:
				return GreedySelectAlignments(alignments, alignmentScoreCompare);
			case GreedyE:
				return GreedySelectAlignments(alignments, std::bind(alignmentECompare, std::placeholders::_1, std::placeholders::_2, options.graphSize, options.readSize));
			case ScheduleInverseESum:
				return ScheduleSelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem aln) { return 1.0 / Evalue(options.graphSize, options.readSize, aln); });
			case ScheduleInverseEProduct:
				return ScheduleSelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem aln) { return -log(Evalue(options.graphSize, options.readSize, aln)); });
			case ScheduleScore:
				return ScheduleSelectAlignments(alignments, S);
			case ScheduleLength:
				return ScheduleSelectAlignments(alignments, [](const AlignmentResult::AlignmentItem aln) { return (aln.alignmentEnd - aln.alignmentStart) + 0.5 - 0.5 / (aln.alignmentScore); });
			default:
			case All:
				return alignments;
		}
		assert(false);
		return alignments;
	}

	std::vector<AlignmentResult::AlignmentItem> SelectECutoff(const std::vector<AlignmentResult::AlignmentItem>& alignments, size_t m, size_t n, double cutoff)
	{
		std::vector<AlignmentResult::AlignmentItem> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			if (Evalue(m, n, alignments[i]) <= cutoff) result.push_back(alignments[i]);
		}
		return result;
	}

}