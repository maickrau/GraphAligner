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
	
	// AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence) const
	// {
	// 	std::vector<typename NodeSlice<WordSlice>::MapItem> nodesliceMap;
	// 	nodesliceMap.resize(params.graph.NodeSize(), {0, 0, 0});
	// 	auto timeStart = std::chrono::system_clock::now();
	// 	assert(params.graph.finalized);
	// 	auto trace = getBacktraceFullStart(sequence, nodesliceMap);
	// 	auto timeEnd = std::chrono::system_clock::now();
	// 	size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	// 	//failed alignment, don't output
	// 	if (std::get<0>(trace) == std::numeric_limits<ScoreType>::max()) return emptyAlignment(time, std::get<2>(trace));
	// 	if (std::get<1>(trace).size() == 0) return emptyAlignment(time, std::get<2>(trace));
	// 	auto result = traceToAlignment(seq_id, sequence, std::get<0>(trace), std::get<1>(trace), std::get<2>(trace));
	// 	result.alignmentStart = std::get<1>(trace)[0].second;
	// 	result.alignmentEnd = std::get<1>(trace).back().second;
	// 	timeEnd = std::chrono::system_clock::now();
	// 	time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	// 	result.elapsedMilliseconds = time;
	// 	return result;
	// }

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

	std::string getSeq(int nodeId, size_t offset, size_t len) const
	{
		std::string result;
		for (size_t i = 0; i < len; i++)
		{
			size_t offsetHere = offset + i;
			size_t nodeIndex = params.graph.GetUnitigNode(nodeId, offsetHere);
			size_t nodeOffset = offsetHere - params.graph.nodeOffset[nodeIndex];
			result += params.graph.NodeSequences(nodeIndex, nodeOffset);
		}
		return result;
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
		assert(seedHit.matchLen > 4);
		assert(getSeq(forwardNodeId, seedHit.nodeOffset, seedHit.matchLen) == sequence.substr(seedHit.seqPos, seedHit.matchLen));
		if (seedHit.seqPos > 0)
		{
			assert(sequence.size() >= seedHit.seqPos + seedHit.matchLen);
			auto backwardPart = CommonUtils::ReverseComplement(sequence.substr(0, seedHit.seqPos + 2));
			auto reversePos = params.graph.GetReversePosition(forwardNodeId, seedHit.nodeOffset + 2);
			assert(reversePos.first == backwardNodeId);
			result.backward = bvAligner.getTraceFromSeed(backwardPart, backwardNodeId, reversePos.second, reusableState);
		}
		if (seedHit.seqPos + seedHit.matchLen < sequence.size())
		{
			auto forwardPart = sequence.substr(seedHit.seqPos + seedHit.matchLen - 2);
			size_t offset = seedHit.nodeOffset + seedHit.matchLen - 2;
			result.forward = bvAligner.getTraceFromSeed(forwardPart, forwardNodeId, offset, reusableState);
			for (auto& item : result.forward.trace)
			{
				item.seqPos += seedHit.seqPos + seedHit.matchLen - 2;
			}
		}
		return result;
	}

	void fixReverseTraceSeqPosAndOrder(std::vector<MatrixPosition>& trace, LengthType end) const
	{
		std::reverse(trace.begin(), trace.end());
		for (size_t i = 0; i < trace.size(); i++)
		{
			trace[i].seqPos = end - trace[i].seqPos;
		}
	}

	AlignmentResult::AlignmentItem getAlignmentFromSeed(const std::string& seq_id, const std::string& sequence, SeedHit seedHit, AlignerGraphsizedState& reusableState) const
	{
		assert(params.graph.finalized);
		assert(seedHit.matchLen >= 4);
		auto timeStart = std::chrono::system_clock::now();

		auto trace = getTwoDirectionalTrace(sequence, seedHit, reusableState);

#ifndef NDEBUG
		if (trace.forward.trace.size() > 0) verifyTrace(trace.forward.trace, sequence, trace.forward.score);
		if (trace.backward.trace.size() > 0) verifyTrace(trace.backward.trace, sequence, trace.backward.score);
#endif
		fixReverseTraceSeqPosAndOrder(trace.backward.trace, seedHit.seqPos + 1);


		//failed alignment, don't output
		if (trace.forward.failed() && trace.backward.failed())
		{
			return VGAlignment::emptyAlignment(0, 0);
		}

		// auto traceVector = getTraceInfo(sequence, trace.backward.trace, trace.forward.trace);

		auto fwresult = VGAlignment::traceToAlignment(params, seq_id, sequence, trace.forward.score, trace.forward.trace, 0, false);
		auto bwresult = VGAlignment::traceToAlignment(params, seq_id, sequence, trace.backward.score, trace.backward.trace, 0, true);
		//failed alignment, don't output
		if (fwresult.alignmentFailed() && bwresult.alignmentFailed())
		{
			return VGAlignment::emptyAlignment(0, 0);
		}
		auto result = VGAlignment::mergeAlignments(params, bwresult, seedHit, fwresult, sequence);
		LengthType seqstart = 0;
		LengthType seqend = 0;
		assert(trace.forward.trace.size() > 0 || trace.backward.trace.size() > 0);
		if (trace.backward.trace.size() > 0 && trace.forward.trace.size() > 0)
		{
			seqstart = trace.backward.trace[0].seqPos;
			seqend = trace.forward.trace.back().seqPos;
		}
		else if (trace.backward.trace.size() > 0)
		{
			seqstart = trace.backward.trace[0].seqPos;
			seqend = trace.backward.trace.back().seqPos;
		}
		else if (trace.forward.trace.size() > 0)
		{
			seqstart = trace.forward.trace[0].seqPos;
			seqend = trace.forward.trace.back().seqPos;
		}
		else
		{
			assert(false);
		}
		result.alignment.set_sequence(sequence.substr(seqstart, seqend - seqstart));
		// result.trace = traceVector;
		result.alignment.set_query_position(seqstart);
		result.alignmentStart = seqstart;
		result.alignmentEnd = seqend;
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
		for (size_t i = 1; i < trace.size(); i++)
		{
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