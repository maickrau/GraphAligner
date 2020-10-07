#include <algorithm>
#include <cstdint>
#include <cmath>
#include "AlignmentSelection.h"
#include "EValue.h"

//an overlap which is larger than the fraction cutoff of the smaller alignment means the alignments are incompatible
//eg alignments 12000bp and 15000bp, overlap of 12000*0.05 = 600bp means they are incompatible
const float OverlapIncompatibleFractionCutoff = 0.05;

namespace AlignmentSelection
{
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
	bool alignmentECompare(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right, size_t m, size_t n, const EValueCalculator& EValueCalc)
	{
		return EValueCalc.getEValue(m, n, left.alignmentLength(), left.alignmentScore) < EValueCalc.getEValue(m, n, right.alignmentLength(), right.alignmentScore);
	}

	bool alignmentScoreCompare(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right, const EValueCalculator& EValueCalc)
	{
		return EValueCalc.getAlignmentScore(left.alignmentLength(), left.alignmentScore) > EValueCalc.getAlignmentScore(right.alignmentLength(), right.alignmentScore);
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
			filteredByE = SelectECutoff(allAlignments, options.graphSize, options.readSize, options.ECutoff, options.EValueCalc);
		}
		std::vector<AlignmentResult::AlignmentItem> filteredByAlignmentScoreFraction;
		if (options.AlignmentScoreFractionCutoff != 0)
		{
			filteredByAlignmentScoreFraction = SelectAlignmentFractionCutoff((options.ECutoff != -1) ? filteredByE : allAlignments, options.AlignmentScoreFractionCutoff, options.EValueCalc);
		}
		const std::vector<AlignmentResult::AlignmentItem>& alignments { options.AlignmentScoreFractionCutoff != 0 ? filteredByAlignmentScoreFraction : ((options.ECutoff != -1) ? filteredByE : allAlignments) };
		switch(options.method)
		{
			case GreedyLength:
				return GreedySelectAlignments(alignments, alignmentLengthCompare);
			case GreedyScore:
				return GreedySelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return alignmentScoreCompare(left, right, options.EValueCalc); });
			case GreedyE:
				return GreedySelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) {return alignmentECompare(left, right, options.graphSize, options.readSize, options.EValueCalc); });
			case ScheduleInverseESum:
				return ScheduleSelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem aln) { return 1.0 / options.EValueCalc.getEValue(options.graphSize, options.readSize, aln.alignmentLength(), aln.alignmentScore); });
			case ScheduleInverseEProduct:
				return ScheduleSelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem aln) { return -log(options.EValueCalc.getEValue(options.graphSize, options.readSize, aln.alignmentLength(), aln.alignmentScore)); });
			case ScheduleScore:
				return ScheduleSelectAlignments(alignments, [options](const AlignmentResult::AlignmentItem aln) { return options.EValueCalc.getAlignmentScore(aln.alignmentLength(), aln.alignmentScore); });
			case ScheduleLength:
				return ScheduleSelectAlignments(alignments, [](const AlignmentResult::AlignmentItem aln) { return (aln.alignmentEnd - aln.alignmentStart) + 0.5 - 0.5 / (aln.alignmentScore); });
			default:
			case All:
				return alignments;
		}
		assert(false);
		return alignments;
	}

	std::vector<AlignmentResult::AlignmentItem> SelectECutoff(const std::vector<AlignmentResult::AlignmentItem>& alignments, size_t m, size_t n, double cutoff, const EValueCalculator& EValueCalc)
	{
		std::vector<AlignmentResult::AlignmentItem> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			if (EValueCalc.getEValue(m, n, alignments[i].alignmentLength(), alignments[i].alignmentScore) <= cutoff) result.push_back(alignments[i]);
		}
		return result;
	}

	std::vector<AlignmentResult::AlignmentItem> SelectAlignmentFractionCutoff(const std::vector<AlignmentResult::AlignmentItem>& alignments, double fraction, const EValueCalculator& EValueCalc)
	{
		std::vector<AlignmentResult::AlignmentItem> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			bool skipped = false;
			for (size_t j = 0; j < alignments.size(); j++)
			{
				if (i == j) continue;
				if (!alignmentIncompatible(alignments[i], alignments[j])) continue;
				if (alignments[i].alignmentXScore < alignments[j].alignmentXScore * fraction)
				{
					skipped = true;
					break;
				}
			}
			if (skipped) continue;
			result.push_back(alignments[i]);
		}
		return result;
	}

	void AddMappingQualities(std::vector<AlignmentResult::AlignmentItem>& alignments)
	{
		for (size_t i = 0; i < alignments.size(); i++)
		{
			assert(alignments[i].alignmentXScore != -1);
			double otherSum = 0;
			for (size_t j = 0; j < alignments.size(); j++)
			{
				if (i == j) continue;
				if (!alignmentIncompatible(alignments[i], alignments[j])) continue;
				assert(alignments[j].alignmentXScore != -1);
				if (alignments[j].alignmentXScore >= alignments[i].alignmentXScore+1)
				{
					otherSum += 10;
					break;
				}
				assert(alignments[j].alignmentXScore <= alignments[i].alignmentXScore+1);
				otherSum += pow(10.0, alignments[j].alignmentXScore - alignments[i].alignmentXScore);
				if (otherSum >= 10) break;
			}
			if (otherSum >= 10)
			{
				alignments[i].mappingQuality = 0;
			}
			else if (otherSum <= 0.000001)
			{
				alignments[i].mappingQuality = 60;
			}
			else
			{
				assert(otherSum >= 0.000001);
				assert(otherSum <= 10);
				alignments[i].mappingQuality = -log(1.0 - 1.0/(1.0 + otherSum)) * 10;
				if (alignments[i].mappingQuality >= 60) alignments[i].mappingQuality = 60;
			}
			assert(alignments[i].mappingQuality >= 0);
			assert(alignments[i].mappingQuality <= 60);
		}
	}


}