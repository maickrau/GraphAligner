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

	std::vector<AlignmentResult::AlignmentItem> SelectAlignments(const std::vector<AlignmentResult::AlignmentItem>& allAlignments, SelectionOptions options)
	{
		// roundabout to fit the signature of const ref while allowing filtering
		std::vector<AlignmentResult::AlignmentItem> filtered;
		bool wasFiltered = false;
		if (options.minAlignmentScore > std::numeric_limits<int>::min())
		{
			filtered = SelectAlignmentScore(wasFiltered ? filtered : allAlignments, options.minAlignmentScore, options.EValueCalc);
			wasFiltered = true;
		}
		if (options.ECutoff != -1)
		{
			filtered = SelectECutoff(wasFiltered ? filtered : allAlignments, options.graphSize, options.readSize, options.ECutoff, options.EValueCalc);
			wasFiltered = true;
		}
		if (options.AlignmentScoreFractionCutoff != 0)
		{
			filtered = SelectAlignmentFractionCutoff(wasFiltered ? filtered : allAlignments, options.AlignmentScoreFractionCutoff, options.EValueCalc);
			wasFiltered = true;
		}
		const std::vector<AlignmentResult::AlignmentItem>& alignments { wasFiltered ? filtered : allAlignments };
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

	std::vector<AlignmentResult::AlignmentItem> SelectAlignmentScore(const std::vector<AlignmentResult::AlignmentItem>& alignments, double cutoff, const EValueCalculator& EValueCalc)
	{
		std::vector<AlignmentResult::AlignmentItem> result;
		for (size_t i = 0; i < alignments.size(); i++)
		{
			if (alignments[i].alignmentXScore >= cutoff) result.push_back(alignments[i]);
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

	std::tuple<size_t, size_t, std::vector<int>> getAlignmentPath(const AlignmentGraph& graph, const std::vector<GraphAlignerCommon<size_t, int32_t, uint64_t>::TraceItem>& trace)
	{
		size_t leftClip = 0;
		size_t rightClip = 0;
		leftClip = trace[0].DPposition.nodeOffset;
		std::vector<int> path;
		path.push_back(trace[0].DPposition.node);
		int currentNodeId = trace[0].DPposition.node;
		size_t currentNodeOffset = trace[0].DPposition.nodeOffset;
		for (size_t pos = 1; pos < trace.size(); pos++)
		{
			int newNodeId = trace[pos].DPposition.node;
			size_t newNodeOffset = trace[pos].DPposition.nodeOffset;
			bool insideNode = !trace[pos-1].nodeSwitch || (newNodeId == currentNodeId && newNodeOffset > currentNodeOffset);
			if (!insideNode)
			{
				currentNodeId = newNodeId;
				currentNodeOffset = newNodeOffset;
				path.push_back(currentNodeId);
			}
		}
		rightClip = graph.OriginalNodeSize(trace.back().DPposition.node) - trace.back().DPposition.nodeOffset;
		return std::make_tuple(leftClip, rightClip, path);
	}

	void RemoveDuplicateAlignments(const AlignmentGraph& graph, std::vector<AlignmentResult::AlignmentItem>& alignments)
	{
		if (alignments.size() <= 1) return;
		std::sort(alignments.begin(), alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });
		std::vector<bool> remove;
		remove.resize(alignments.size(), false);
		size_t blockStart = 0;
		for (size_t i = 1; i <= alignments.size(); i++)
		{
			if (i < alignments.size() && alignments[i].alignmentStart == alignments[i-1].alignmentStart && alignments[i].alignmentEnd == alignments[i-1].alignmentEnd) continue;
			if (i == blockStart + 1)
			{
				blockStart = i;
				continue;
			}
			std::vector<std::tuple<size_t, size_t, std::vector<int>>> paths;
			std::vector<size_t> pathOrder;
			for (size_t j = blockStart; j < i; j++)
			{
				paths.emplace_back(getAlignmentPath(graph, alignments[j].trace->trace));
				pathOrder.push_back(j-blockStart);
			}
			std::sort(pathOrder.begin(), pathOrder.end(), [&paths](size_t left, size_t right) { return paths[left] < paths[right]; });
			size_t pathOrderBlockStart = 0;
			for (size_t j = 1; j <= pathOrder.size(); j++)
			{
				if (j < pathOrder.size() && paths[j-1] == paths[j]) continue;
				size_t maxAlignmentScoreIndex = pathOrderBlockStart;
				for (size_t k = pathOrderBlockStart+1; k < j; k++)
				{
					if (alignments[blockStart+pathOrder[k]].alignmentXScore > alignments[blockStart+pathOrder[maxAlignmentScoreIndex]].alignmentXScore)
					{
						maxAlignmentScoreIndex = k;
					}
				}
				for (size_t k = pathOrderBlockStart; k < j; k++)
				{
					remove[blockStart+pathOrder[k]] = true;
				}
				remove[blockStart+pathOrder[maxAlignmentScoreIndex]] = false;
				pathOrderBlockStart = j;
			}
			blockStart = i;
		}
		for (size_t i = alignments.size()-1; i < alignments.size(); i--)
		{
			if (remove[i])
			{
				if (i != alignments.size()-1) std::swap(alignments[i], alignments.back());
				alignments.pop_back();
			}
		}
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