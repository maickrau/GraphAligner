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
	enum SelectionMethod
	{
		GreedyLength,
		GreedyScore,
		GreedyE,
		ScheduleInverseESum,
		ScheduleInverseEProduct,
		ScheduleScore,
		ScheduleLength,
		All
	};
	struct SelectionOptions
	{
		SelectionMethod method;
		size_t graphSize;
		size_t readSize;
		double ECutoff;
		EValueCalculator EValueCalc;
	};
	std::vector<AlignmentResult::AlignmentItem> SelectAlignments(const std::vector<AlignmentResult::AlignmentItem>& alignments, SelectionOptions options);
	bool alignmentIncompatible(const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right);

	template <typename AlnScorer>
	std::vector<AlignmentResult::AlignmentItem> GreedySelectAlignments(const std::vector<AlignmentResult::AlignmentItem>& alignments, AlnScorer alnScorer)
	{
		std::vector<size_t> items;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			items.push_back(i);
		}
		std::sort(items.begin(), items.end(), [&alignments, alnScorer](size_t left, size_t right) { return alnScorer(alignments[left], alignments[right]); });
		std::vector<AlignmentResult::AlignmentItem> result;
		for (auto i : items)
		{
			if (!std::any_of(result.begin(), result.end(), [&alignments, i](const AlignmentResult::AlignmentItem& existing) { return alignmentIncompatible(existing, alignments[i]); }))
			{
				result.push_back(alignments[i]);
			}
		}
		return result;
	}

	template <typename AlnScorer>
	std::vector<AlignmentResult::AlignmentItem> ScheduleSelectAlignments(const std::vector<AlignmentResult::AlignmentItem>& alignments, AlnScorer alnScorer)
	{
		std::vector<size_t> items;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			items.push_back(i);
		}
		std::sort(items.begin(), items.end(), [&alignments](size_t left, size_t right) { return alignments[left].alignmentEnd < alignments[right].alignmentEnd; });
		std::vector<size_t> backtrace;
		std::vector<double> score;
		backtrace.resize(items.size(), std::numeric_limits<size_t>::max());
		score.resize(items.size(), 0);
		for (size_t i = 0; i < items.size(); i++)
		{
			double rawScore = alnScorer(alignments[items[i]]);
			score[i] = rawScore;
			for (size_t j = 0; j < i; j++)
			{
				if (alignmentIncompatible(alignments[items[i]], alignments[items[j]])) continue;
				if (score[j] + rawScore > score[i])
				{
					backtrace[i] = j;
					score[i] = score[j] + rawScore;
				}
			}
		}
		size_t maxPos = 0;
		for (size_t i = 0; i < items.size(); i++)
		{
			if (score[i] > score[maxPos]) maxPos = i;
		}
		std::vector<AlignmentResult::AlignmentItem> result;
		while (maxPos != std::numeric_limits<size_t>::max())
		{
			result.push_back(alignments[items[maxPos]]);
			assert(backtrace[maxPos] < maxPos || backtrace[maxPos] == std::numeric_limits<size_t>::max());
			maxPos = backtrace[maxPos];
		}
		return result;
	}

	std::vector<AlignmentResult::AlignmentItem> SelectECutoff(const std::vector<AlignmentResult::AlignmentItem>& alignments, size_t m, size_t n, double cutoff, const EValueCalculator& EValueCalc);

};

#endif
