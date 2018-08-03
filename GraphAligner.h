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
	using SeedHit = typename Common::SeedHit;
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
		alnItem.alignmentStart = trace.trace[0].seqPos;
		alnItem.alignmentEnd = trace.trace.back().seqPos;
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
			std::string seedInfo = std::to_string(seedHits[i].nodeID) + (seedHits[i].reverse ? "-" : "+") + "," + std::to_string(seedHits[i].seqPos);
			logger << seq_id << " seed " << i << "/" << seedHits.size() << " " << seedInfo;
			assertSetRead(seq_id, seedInfo);
			// auto nodeIndex = params.graph.nodeLookup.at(std::get<0>(seedHits[i]) * 2);
			// auto pos = std::get<1>(seedHits[i]);
			// if (std::any_of(triedAlignmentNodes.begin(), triedAlignmentNodes.end(), [nodeIndex, pos](auto triple) { return std::get<0>(triple) <= pos && std::get<1>(triple) >= pos && std::get<2>(triple) == nodeIndex; }))
			// {
			// 	logger << "seed " << i << " already aligned" << BufferedWriter::Flush;
			// 	continue;
			// }
			logger << BufferedWriter::Flush;
			auto item = getAlignmentFromSeed(seq_id, sequence, seedHits[i], reusableState);
			if (item.alignmentFailed()) continue;
			result.alignments.push_back(item);
			// addAlignmentNodes(triedAlignmentNodes, item);
			if (params.sloppyOptimizations && item.alignmentStart == 0 && item.alignmentEnd >= sequence.size() - params.graph.DBGOverlap - 1)
			{
				break;
			}
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
			assert(result.backward.trace.back().seqPos == -1 && params.graph.nodeIDs[result.backward.trace.back().node] == backwardNodeId && params.graph.nodeOffset[result.backward.trace.back().node] + result.backward.trace.back().nodeOffset == reversePos.second);
			result.backward.trace.pop_back();
			std::reverse(result.backward.trace.begin(), result.backward.trace.end());
		}
		if (!result.forward.failed())
		{
			assert(result.forward.trace.back().seqPos == -1 && params.graph.nodeIDs[result.forward.trace.back().node] == forwardNodeId && params.graph.nodeOffset[result.forward.trace.back().node] + result.forward.trace.back().nodeOffset == seedHit.nodeOffset);
			result.forward.trace.pop_back();
			std::reverse(result.forward.trace.begin(), result.forward.trace.end());
		}
		return result;
	}

	void fixForwardTraceSeqPos(std::vector<MatrixPosition>& trace, LengthType start) const
	{
		for (size_t i = 0; i < trace.size(); i++)
		{
			trace[i].seqPos += start;
		}
	}

	void fixReverseTraceSeqPosAndOrder(std::vector<MatrixPosition>& trace, LengthType end) const
	{
		std::reverse(trace.begin(), trace.end());
		for (size_t i = 0; i < trace.size(); i++)
		{
			assert(trace[i].seqPos <= end || trace[i].seqPos == -1);
			trace[i].seqPos = end - trace[i].seqPos;
			size_t offset = params.graph.nodeOffset[trace[i].node] + trace[i].nodeOffset;
			auto reversePos = params.graph.GetReversePosition(params.graph.nodeIDs[trace[i].node], offset);
			trace[i].node = params.graph.GetUnitigNode(reversePos.first, reversePos.second);
			trace[i].nodeOffset = reversePos.second - params.graph.nodeOffset[trace[i].node];
			assert(trace[i].nodeOffset < params.graph.NodeLength(trace[i].node));
		}
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

// #ifndef NDEBUG
// 		if (!trace.forward.failed() && !trace.backward.failed())
// 		{
// 			assert(trace.backward.trace.back().seqPos == trace.forward.trace[0].seqPos - 1);
// 			auto debugBwPos = trace.backward.trace.back();
// 			auto debugBwNodeId = params.graph.nodeIDs[debugBwPos.node];
// 			auto debugBwOffset = params.graph.nodeOffset[debugBwPos.node] + debugBwPos.nodeOffset;
// 			auto debugOldNodeidPos = params.graph.GetReversePosition(debugBwNodeId, debugBwOffset);
// 			auto debugOldNode = params.graph.GetUnitigNode(debugOldNodeidPos.first, debugOldNodeidPos.second);
// 			auto debugOldOffset = debugOldNodeidPos.second - params.graph.nodeOffset[debugOldNode];
// 			assert(debugOldOffset >= 0);
// 			assert(debugOldOffset < params.graph.NodeLength(debugOldNode));
// 			if (trace.forward.trace[0].nodeOffset > 0)
// 			{
// 				assert(debugOldNode == trace.forward.trace[0].node);
// 				assert(debugOldOffset <= trace.forward.trace[0].nodeOffset);
// 			}
// 			else
// 			{
// 				bool foundNeighbor = false;
// 				for (auto neighbor : params.graph.inNeighbors[trace.forward.trace[0].node])
// 				{
// 					if (debugOldNode == neighbor && debugOldOffset == params.graph.NodeLength(neighbor)-1) foundNeighbor = true;
// 				}
// 				assert(foundNeighbor || (debugOldNode == trace.forward.trace[0].node && debugOldOffset == 0));
// 			}
// 		}
// #endif

		//failed alignment, don't output
		if (trace.forward.failed() && trace.backward.failed())
		{
			return VGAlignment::emptyAlignment(0, 0);
		}

		// auto traceVector = getTraceInfo(sequence, trace.backward.trace, trace.forward.trace);

		auto mergedTrace = trace.backward;
		if (trace.backward.failed()) mergedTrace.score = 0;

		int seedHitNodeId;
		if (seedHit.reverse)
		{
			seedHitNodeId = seedHit.nodeID * 2 + 1;
		}
		else
		{
			seedHitNodeId = seedHit.nodeID * 2;
		}
		auto middleNode = params.graph.GetUnitigNode(seedHitNodeId, seedHit.nodeOffset);
		mergedTrace.trace.emplace_back(middleNode, seedHit.nodeOffset - params.graph.nodeOffset[middleNode], seedHit.seqPos);
		mergedTrace.trace.insert(mergedTrace.trace.end(), trace.forward.trace.begin(), trace.forward.trace.end());
		if (!trace.forward.failed()) mergedTrace.score += trace.forward.score;

		auto result = VGAlignment::traceToAlignment(params, seq_id, sequence, mergedTrace.score, mergedTrace.trace, 0, false);

		assert(!result.alignmentFailed());
		LengthType seqstart = 0;
		LengthType seqend = 0;
		assert(mergedTrace.trace.size() > 0);
		seqstart = mergedTrace.trace[0].seqPos;
		seqend = mergedTrace.trace.back().seqPos;
		assert(seqend < sequence.size());
		result.alignment.set_sequence(sequence.substr(seqstart, seqend - seqstart + 1));
		// result.trace = traceVector;
		result.alignment.set_query_position(seqstart);
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
	void verifyTrace(const std::vector<MatrixPosition>& trace, const std::string& sequence, volatile ScoreType score) const
	{
		size_t start = 0;
		while (trace[start].seqPos == -1)
		{
			start++;
			assert(start < trace.size());
		}
		start++;
		for (size_t i = start; i < trace.size(); i++)
		{
			assert(trace[i].seqPos < sequence.size());
			auto newpos = trace[i];
			auto oldpos = trace[i-1];
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