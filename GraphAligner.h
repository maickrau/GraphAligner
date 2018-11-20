#ifndef GraphAligner_H
#define GraphAligner_H

#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include "AlignmentGraph.h"
#include "CommonUtils.h"
#include "GraphAlignerWrapper.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerCommon.h"
#include "GraphAlignerVGAlignment.h"
#include "GraphAlignerBitvectorBanded.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAligner
{
private:
	using VGAlignment = GraphAlignerVGAlignment<LengthType, ScoreType, Word>;
	using BitvectorAligner = GraphAlignerBitvectorBanded<LengthType, ScoreType, Word>;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	BitvectorAligner bvAligner;
	mutable BufferedWriter logger;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	const Params& params;
public:

	GraphAligner(const Params& params) :
	bvAligner(params),
	logger(),
	params(params)
	{
		if (!params.quietMode) logger = { std::cerr };
	}
	
	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, AlignerGraphsizedState& reusableState) const
	{
		AlignmentResult result;
		auto timeStart = std::chrono::system_clock::now();
		assert(params.graph.finalized);
		auto trace = getBacktraceFullStart(sequence, reusableState);
		auto timeEnd = std::chrono::system_clock::now();
		size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		//failed alignment, don't output
		if (trace.score == std::numeric_limits<ScoreType>::max()) return result;
		if (trace.trace.size() == 0) return result;
		auto alnItem = VGAlignment::traceToAlignment(params, seq_id, sequence, trace.score, trace.trace, 0, false);
		alnItem.alignmentStart = trace.trace[0].first.seqPos;
		alnItem.alignmentEnd = trace.trace.back().first.seqPos;
		timeEnd = std::chrono::system_clock::now();
		time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		alnItem.elapsedMilliseconds = time;
		result.alignments.push_back(alnItem);
		return result;
	}

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, const std::vector<SeedHit>& seedHits, AlignerGraphsizedState& reusableState) const
	{
		assert(params.graph.finalized);
		AlignmentResult result;
		assert(seedHits.size() > 0);
		// std::vector<std::tuple<size_t, size_t, size_t>> triedAlignmentNodes;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			std::string seedInfo = std::to_string(seedHits[i].nodeID) + (seedHits[i].reverse ? "-" : "+") + "," + std::to_string(seedHits[i].seqPos) + "," + std::to_string(seedHits[i].matchLen) + "," + std::to_string(seedHits[i].nodeOffset);
			logger << seq_id << " seed " << i << "/" << seedHits.size() << " " << seedInfo;
			assertSetRead(seq_id, seedInfo);
			if (params.sloppyOptimizations)
			{
				bool found = false;
				for (auto aln : result.alignments)
				{
					if (aln.alignmentStart <= seedHits[i].seqPos && aln.alignmentEnd >= seedHits[i].seqPos)
					{
						logger << " skipped";
						logger << BufferedWriter::Flush;
						found = true;
						break;
					}
				}
				if (found) continue;
			}
			logger << BufferedWriter::Flush;
			auto item = getAlignmentFromSeed(seq_id, sequence, seedHits[i], reusableState);
			if (item.alignmentFailed()) continue;
			result.alignments.push_back(item);
		}
		assertSetRead(seq_id, "No seed");

		return result;
	}

private:

	OnewayTrace getBacktraceFullStart(const std::string& sequence, AlignerGraphsizedState& reusableState) const
	{
		return bvAligner.getBacktraceFullStart(sequence, reusableState);
	}

	Trace getTwoDirectionalTrace(const std::string& sequence, SeedHit seedHit, AlignerGraphsizedState& reusableState) const
	{
		assert(seedHit.seqPos >= 0);
		assert(seedHit.seqPos < sequence.size());
		int forwardNodeId;
		int backwardNodeId;
		if (seedHit.reverse)
		{
			forwardNodeId = seedHit.nodeID * 2 + 1;
			backwardNodeId = seedHit.nodeID * 2;
		}
		else
		{
			forwardNodeId = seedHit.nodeID * 2;
			backwardNodeId = seedHit.nodeID * 2 + 1;
		}
		Trace result;
		result.backward.score = std::numeric_limits<ScoreType>::max();
		result.forward.score = std::numeric_limits<ScoreType>::max();
		if (seedHit.seqPos > 0)
		{
			auto backwardPart = CommonUtils::ReverseComplement(sequence.substr(0, seedHit.seqPos));
			auto reversePos = params.graph.GetReversePosition(forwardNodeId, seedHit.nodeOffset);
			assert(reversePos.first == backwardNodeId);
			result.backward = bvAligner.getReverseTraceFromSeed(backwardPart, backwardNodeId, reversePos.second, reusableState);
		}
		if (seedHit.seqPos < sequence.size()-1)
		{
			auto forwardPart = sequence.substr(seedHit.seqPos+1);
			size_t offset = seedHit.nodeOffset;
			result.forward = bvAligner.getReverseTraceFromSeed(forwardPart, forwardNodeId, offset, reusableState);
		}

		if (!result.backward.failed())
		{
			auto reversePos = params.graph.GetReversePosition(forwardNodeId, seedHit.nodeOffset);
			assert(result.backward.trace.back().first.seqPos == (size_t)-1 && params.graph.nodeIDs[result.backward.trace.back().first.node] == backwardNodeId && params.graph.nodeOffset[result.backward.trace.back().first.node] + result.backward.trace.back().first.nodeOffset == reversePos.second);
			std::reverse(result.backward.trace.begin(), result.backward.trace.end());
		}
		if (!result.forward.failed())
		{
			assert(result.forward.trace.back().first.seqPos == (size_t)-1 && params.graph.nodeIDs[result.forward.trace.back().first.node] == forwardNodeId && params.graph.nodeOffset[result.forward.trace.back().first.node] + result.forward.trace.back().first.nodeOffset == seedHit.nodeOffset);
			std::reverse(result.forward.trace.begin(), result.forward.trace.end());
		}
		return result;
	}

	void fixForwardTraceSeqPos(std::vector<std::pair<MatrixPosition, bool>>& trace, LengthType start) const
	{
		for (size_t i = 0; i < trace.size(); i++)
		{
			trace[i].first.seqPos += start;
			auto nodeIndex = trace[i].first.node;
			trace[i].first.node = params.graph.nodeIDs[nodeIndex];
			trace[i].first.nodeOffset += params.graph.nodeOffset[nodeIndex];
		}
	}

	void fixReverseTraceSeqPosAndOrder(std::vector<std::pair<MatrixPosition, bool>>& trace, LengthType end) const
	{
		if (trace.size() == 0) return;
		std::reverse(trace.begin(), trace.end());
		for (size_t i = 0; i < trace.size(); i++)
		{
			assert(trace[i].first.seqPos <= end || trace[i].first.seqPos == (size_t)-1);
			trace[i].first.seqPos = end - trace[i].first.seqPos;
			size_t offset = params.graph.nodeOffset[trace[i].first.node] + trace[i].first.nodeOffset;
			auto reversePos = params.graph.GetReversePosition(params.graph.nodeIDs[trace[i].first.node], offset);
			assert(reversePos.second < params.graph.originalNodeSize.at(params.graph.nodeIDs[trace[i].first.node]));
			trace[i].first.node = reversePos.first;
			trace[i].first.nodeOffset = reversePos.second;
		}
		for (size_t i = 0; i < trace.size() - 1; i++)
		{
			trace[i].second = trace[i+1].second;
		}
		trace.back().second = false;
	}

	AlignmentResult::AlignmentItem getAlignmentFromSeed(const std::string& seq_id, const std::string& sequence, SeedHit seedHit, AlignerGraphsizedState& reusableState) const
	{
		assert(params.graph.finalized);
		auto timeStart = std::chrono::system_clock::now();

		auto trace = getTwoDirectionalTrace(sequence, seedHit, reusableState);

#ifndef NDEBUG
		if (trace.forward.trace.size() > 0) verifyTrace(trace.forward.trace, sequence, trace.forward.score);
		if (trace.backward.trace.size() > 0) verifyTrace(trace.backward.trace, sequence, trace.backward.score);
#endif
		fixReverseTraceSeqPosAndOrder(trace.backward.trace, seedHit.seqPos-1);
		fixForwardTraceSeqPos(trace.forward.trace, seedHit.seqPos+1);

		//failed alignment, don't output
		if (trace.forward.failed() && trace.backward.failed())
		{
			return VGAlignment::emptyAlignment(0, 0);
		}

		// auto traceVector = getTraceInfo(sequence, trace.backward.trace, trace.forward.trace);

		auto mergedTrace = trace.backward;
		if (trace.backward.failed()) mergedTrace.score = 0;

		if (!trace.backward.failed() && !trace.forward.failed())
		{
			assert(mergedTrace.trace.back().first == trace.forward.trace[0].first);
			mergedTrace.trace.pop_back();
			mergedTrace.trace.insert(mergedTrace.trace.end(), trace.forward.trace.begin(), trace.forward.trace.end());
			mergedTrace.score += trace.forward.score;
		}
		else if (!trace.forward.failed())
		{
			mergedTrace.trace.insert(mergedTrace.trace.end(), trace.forward.trace.begin(), trace.forward.trace.end());
			mergedTrace.score += trace.forward.score;
		}

		auto result = VGAlignment::traceToAlignment(params, seq_id, sequence, mergedTrace.score, mergedTrace.trace, 0, false);

		assert(!result.alignmentFailed());
		LengthType seqstart = 0;
		LengthType seqend = 0;
		assert(mergedTrace.trace.size() > 0);
		seqstart = mergedTrace.trace[0].first.seqPos;
		seqend = mergedTrace.trace.back().first.seqPos;
		assert(seqend < sequence.size());
		result.alignment->set_sequence(sequence.substr(seqstart, seqend - seqstart + 1));
		// result.trace = traceVector;
		result.alignment->set_query_position(seqstart);
		result.alignmentStart = seqstart;
		result.alignmentEnd = seqend + 1;
		auto timeEnd = std::chrono::system_clock::now();
		size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		result.elapsedMilliseconds = time;
		return result;
	}

	// void addAlignmentNodes(std::vector<std::tuple<size_t, size_t, size_t>>& tried, const AlignmentResult::AlignmentItem& trace) const
	// {
	// 	assert(trace.trace.size() > 0);
	// 	LengthType currentNode = trace.trace[0].nodeID;
	// 	size_t currentReadStart = trace.trace[0].readpos;
	// 	for (size_t i = 1; i < trace.trace.size(); i++)
	// 	{
	// 		if (trace.trace[i].nodeID != currentNode)
	// 		{
	// 			tried.emplace_back(currentReadStart, trace.trace[i-1].readpos, currentNode);
	// 			currentNode = trace.trace[i].nodeID;
	// 			currentReadStart = trace.trace[i].readpos;
	// 		}
	// 	}
	// }

#ifndef NDEBUG
	void verifyTrace(const std::vector<std::pair<MatrixPosition, bool>>& trace, const std::string& sequence, volatile ScoreType score) const
	{
		size_t start = 0;
		while (trace[start].first.seqPos == (size_t)-1)
		{
			start++;
			assert(start < trace.size());
		}
		start++;
		for (size_t i = start; i < trace.size(); i++)
		{
			assert(trace[i].first.seqPos < sequence.size());
			auto newpos = trace[i].first;
			auto oldpos = trace[i-1].first;
			auto oldNodeIndex = oldpos.node;
			auto newNodeIndex = newpos.node;
			assert(newpos.seqPos != oldpos.seqPos || newpos.node != oldpos.node || newpos.nodeOffset != oldpos.nodeOffset);
			if (oldNodeIndex == newNodeIndex)
			{
				assert((newpos.nodeOffset >= oldpos.nodeOffset && newpos.seqPos >= oldpos.seqPos) || (oldpos.nodeOffset == params.graph.NodeLength(newNodeIndex)-1 && newpos.nodeOffset == 0 && (newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1)));
				continue;
			}
			assert(oldNodeIndex != newNodeIndex);
			assert(std::find(params.graph.outNeighbors[oldpos.node].begin(), params.graph.outNeighbors[oldpos.node].end(), newpos.node) != params.graph.outNeighbors[oldpos.node].end());
			assert(newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1);
			assert(oldpos.nodeOffset == params.graph.NodeLength(oldNodeIndex)-1);
			assert(newpos.nodeOffset == 0);
		}
	}
#endif

};

#endif