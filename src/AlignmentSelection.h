#ifndef AlignmentSelection_h
#define AlignmentSelection_h

#include <cmath>
#include <cstdint>
#include <vector>
#include "vg.pb.h"

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
		All,
		ECutoff
	};
	struct SelectionOptions
	{
		SelectionMethod method;
		size_t graphSize;
		size_t readSize;
		double ECutoff;
	};
	std::vector<vg::Alignment> SelectAlignments(std::vector<vg::Alignment> alns, SelectionOptions method);
	std::vector<vg::Alignment*> SelectAlignments(std::vector<vg::Alignment*> alns, SelectionOptions method);
	double S(const vg::Alignment* const aln);
	double Evalue(size_t m, size_t n, const vg::Alignment* const aln);
	bool alignmentIncompatible(const vg::Alignment* const left, const vg::Alignment* const right);
	bool alignmentLengthCompare(const vg::Alignment* const left, const vg::Alignment* const right);
	bool alignmentScoreCompare(const vg::Alignment* const left, const vg::Alignment* const right);
	bool alignmentECompare(const vg::Alignment* const left, const vg::Alignment* const right, size_t m, size_t n);

	template <typename T, typename AlnGetter, typename AlnScorer>
	std::vector<T> GreedySelectAlignments(std::vector<T> alignments, AlnGetter alnGetter, AlnScorer alnScorer)
	{
		std::function<const vg::Alignment*(const T&)> f = [alnGetter](const T& aln) { return (const vg::Alignment*)alnGetter(aln); };
		std::sort(alignments.begin(), alignments.end(), [f, alnScorer](const T& left, const T& right) { return alnScorer(f(left), f(right)); });
		std::vector<T> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			const vg::Alignment* const aln = f(alignments[i]);
			if (!std::any_of(result.begin(), result.end(), [aln, f](const T& existing) { return alignmentIncompatible(f(existing), aln); }))
			{
				result.push_back(alignments[i]);
			}
		}
		return result;
	}

	template <typename T, typename AlnGetter, typename AlnScorer>
	std::vector<T> ScheduleSelectAlignments(const std::vector<T>& alignments, AlnGetter alnGetter, AlnScorer alnScorer)
	{
		std::function<const vg::Alignment*(const T&)> f = [alnGetter](const T& aln) { return (const vg::Alignment*)alnGetter(aln); };
		std::vector<std::pair<const vg::Alignment*, size_t>> items;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			const vg::Alignment* aln = f(alignments[i]);
			items.emplace_back(aln, i);
		}
		std::sort(items.begin(), items.end(), [](const std::pair<const vg::Alignment*, size_t>& left, const std::pair<const vg::Alignment*, size_t>& right) { return left.first->query_position() + left.first->sequence().size() < right.first->query_position() + right.first->sequence().size(); });
		std::vector<size_t> backtrace;
		std::vector<double> score;
		backtrace.resize(items.size(), std::numeric_limits<size_t>::max());
		score.resize(items.size(), 0);
		for (size_t i = 0; i < items.size(); i++)
		{
			double rawScore = alnScorer(items[i].first);
			score[i] = rawScore;
			for (size_t j = 0; j < i; j++)
			{
				if (alignmentIncompatible(items[i].first, items[j].first)) continue;
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
		std::vector<T> result;
		while (maxPos != std::numeric_limits<size_t>::max())
		{
			result.push_back(alignments[items[maxPos].second]);
			maxPos = backtrace[maxPos];
		}
		return result;
	}

	template <typename T, typename AlnGetter>
	std::vector<T> SelectECutoff(std::vector<T> alignments, AlnGetter alnGetter, size_t m, size_t n, double cutoff)
	{
		std::function<const vg::Alignment*(const T&)> f = [alnGetter](const T& aln) { return (const vg::Alignment*)alnGetter(aln); };
		std::vector<T> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			const vg::Alignment* const aln = f(alignments[i]);
			if (Evalue(m, n, aln) <= cutoff) result.push_back(alignments[i]);
		}
		return result;
	}

	template <typename T, typename AlnGetter>
	std::vector<T> SelectAlignments(std::vector<T> alignments, SelectionOptions options, AlnGetter alnGetter)
	{
		switch(options.method)
		{
			case GreedyLength:
				return GreedySelectAlignments(alignments, alnGetter, alignmentLengthCompare);
			case GreedyScore:
				return GreedySelectAlignments(alignments, alnGetter, alignmentScoreCompare);
			case GreedyE:
				return GreedySelectAlignments(alignments, alnGetter, std::bind(alignmentECompare, std::placeholders::_1, std::placeholders::_2, options.graphSize, options.readSize));
			case ECutoff:
				return SelectECutoff(alignments, alnGetter, options.graphSize, options.readSize, options.ECutoff);
			case ScheduleInverseESum:
				return ScheduleSelectAlignments(alignments, alnGetter, [options](const vg::Alignment* const aln) { return 1.0 / Evalue(options.graphSize, options.readSize, aln); });
			case ScheduleInverseEProduct:
				return ScheduleSelectAlignments(alignments, alnGetter, [options](const vg::Alignment* const aln) { return -log(Evalue(options.graphSize, options.readSize, aln)); });
			case ScheduleScore:
				return ScheduleSelectAlignments(alignments, alnGetter, S);
			case ScheduleLength:
				return ScheduleSelectAlignments(alignments, alnGetter, [](const vg::Alignment* const aln) { return aln->sequence().size() + 0.5 - 0.5 / aln->score(); });
			default:
			case All:
				return alignments;
		}
		assert(false);
		return alignments;
	}

};

#endif
