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
#include "GraphAlignerBitvector.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAligner
{
private:
	using VGAlignment = GraphAlignerVGAlignment<LengthType, ScoreType, Word>;
	using BitvectorAligner = GraphAlignerBitvector<LengthType, ScoreType, Word>;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using Params = typename Common::Params;
	using SeedHit = typename Common::SeedHit;
	using MatrixPosition = typename Common::MatrixPosition;
	mutable BitvectorAligner bvAligner;
	mutable BufferedWriter logger;
	const Params& params;
public:

	GraphAligner(const Params& params) :
	bvAligner(params),
	logger(std::cerr),
	params(params)
	{
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

	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, const std::vector<SeedHit>& seedHits) const
	{
		assert(params.graph.finalized);
		AlignmentResult result;
		assert(seedHits.size() > 0);
		// std::vector<std::tuple<size_t, size_t, size_t>> triedAlignmentNodes;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			logger << "seed " << i << "/" << seedHits.size() << " " << seedHits[i].nodeID << (seedHits[i].reverse ? "-" : "+") << "," << seedHits[i].seqPos;
			// auto nodeIndex = params.graph.nodeLookup.at(std::get<0>(seedHits[i]) * 2);
			// auto pos = std::get<1>(seedHits[i]);
			// if (std::any_of(triedAlignmentNodes.begin(), triedAlignmentNodes.end(), [nodeIndex, pos](auto triple) { return std::get<0>(triple) <= pos && std::get<1>(triple) >= pos && std::get<2>(triple) == nodeIndex; }))
			// {
			// 	logger << "seed " << i << " already aligned" << BufferedWriter::Flush;
			// 	continue;
			// }
			logger << BufferedWriter::Flush;
			auto item = getAlignmentFromSeed(seq_id, sequence, seedHits[i]);
			if (item.alignmentFailed()) continue;
			// addAlignmentNodes(triedAlignmentNodes, item);
			result.alignments.push_back(item);
		}

		return result;
	}

private:

	AlignmentResult::AlignmentItem getAlignmentFromSeed(const std::string& seq_id, const std::string& sequence, SeedHit seedHit) const
	{
		assert(params.graph.finalized);
		auto timeStart = std::chrono::system_clock::now();

		auto trace = bvAligner.getTraceFromSeed(sequence, seedHit);

#ifndef NDEBUG
		if (trace.forward.trace.size() > 0) verifyTrace(trace.forward.trace, sequence, trace.forward.score);
		if (trace.backward.trace.size() > 0) verifyTrace(trace.backward.trace, sequence, trace.backward.score);
#endif

		//failed alignment, don't output
		if (trace.forward.failed() && trace.backward.failed())
		{
			return VGAlignment::emptyAlignment(0, 0);
		}

		// auto traceVector = getTraceInfo(sequence, trace.backward.trace, trace.forward.trace);

		auto fwresult = VGAlignment::traceToAlignment(params, seq_id, sequence, trace.forward.score, trace.forward.trace, 0);
		auto bwresult = VGAlignment::traceToAlignment(params, seq_id, sequence, trace.backward.score, trace.backward.trace, 0);
		//failed alignment, don't output
		if (fwresult.alignmentFailed() && bwresult.alignmentFailed())
		{
			return VGAlignment::emptyAlignment(0, 0);
		}
		auto result = VGAlignment::mergeAlignments(params, bwresult, fwresult);
		LengthType seqstart = 0;
		LengthType seqend = 0;
		assert(trace.forward.trace.size() > 0 || trace.backward.trace.size() > 0);
		if (trace.backward.trace.size() > 0 && trace.forward.trace.size() > 0)
		{
			seqstart = trace.backward.trace[0].second;
			seqend = trace.forward.trace.back().second;
		}
		else if (trace.backward.trace.size() > 0)
		{
			seqstart = trace.backward.trace[0].second;
			seqend = trace.backward.trace.back().second;
		}
		else if (trace.forward.trace.size() > 0)
		{
			seqstart = trace.forward.trace[0].second;
			seqend = trace.forward.trace.back().second;
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

	std::vector<AlignmentResult::TraceItem> getTraceInfo(const std::string& sequence, const std::vector<MatrixPosition>& bwtrace, const std::vector<MatrixPosition>& fwtrace) const
	{
		std::vector<AlignmentResult::TraceItem> result;
		if (bwtrace.size() > 0)
		{
			auto bw = getTraceInfoInner(sequence, bwtrace);
			result.insert(result.end(), bw.begin(), bw.end());
		}
		if (bwtrace.size() > 0 && fwtrace.size() > 0)
		{
			auto nodeid = params.graph.IndexToNode(fwtrace[0].first);
			result.emplace_back();
			result.back().type = AlignmentResult::TraceMatchType::FORWARDBACKWARDSPLIT;
			result.back().nodeID = params.graph.nodeIDs[nodeid] / 2;
			result.back().reverse = nodeid % 2 == 1;
			result.back().offset = fwtrace[0].first - params.graph.NodeStart(nodeid);
			result.back().readpos = fwtrace[0].second;
			result.back().graphChar = params.graph.NodeSequences(fwtrace[0].first);
			result.back().readChar = sequence[fwtrace[0].second];
		}
		if (fwtrace.size() > 0)
		{
			auto fw = getTraceInfoInner(sequence, fwtrace);
			result.insert(result.end(), fw.begin(), fw.end());
		}
		return result;
	}

	std::vector<AlignmentResult::TraceItem> getTraceInfoInner(const std::string& sequence, const std::vector<MatrixPosition>& trace) const
	{
		std::vector<AlignmentResult::TraceItem> result;
		for (size_t i = 1; i < trace.size(); i++)
		{
			auto newpos = trace[i];
			auto oldpos = trace[i-1];
			assert(newpos.second == oldpos.second || newpos.second == oldpos.second+1);
			assert(newpos.second != oldpos.second || newpos.first != oldpos.first);
			auto oldNodeIndex = params.graph.IndexToNode(oldpos.first);
			auto newNodeIndex = params.graph.IndexToNode(newpos.first);
			if (oldpos.first == params.graph.NodeEnd(oldNodeIndex)-1)
			{
				assert(newpos.first == oldpos.first || newpos.first == params.graph.NodeStart(newNodeIndex));
			}
			else
			{
				assert(newpos.first == oldpos.first || newpos.first == oldpos.first+1);
			}
			bool diagonal = true;
			if (newpos.second == oldpos.second) diagonal = false;
			if (newpos.first == oldpos.first)
			{
				auto newNodeIndex = params.graph.IndexToNode(newpos.first);
				if (newpos.second == oldpos.second+1 && params.graph.NodeEnd(newNodeIndex) == params.graph.NodeStart(newNodeIndex)+1 && std::find(params.graph.outNeighbors[newNodeIndex].begin(), params.graph.outNeighbors[newNodeIndex].end(), newNodeIndex) != params.graph.outNeighbors[newNodeIndex].end())
				{
					//one node self-loop, diagonal is valid
				}
				else
				{
					diagonal = false;
				}
			}
			result.emplace_back();
			result.back().nodeID = params.graph.nodeIDs[newNodeIndex] / 2;
			result.back().reverse = params.graph.nodeIDs[newNodeIndex] % 2 == 1;
			result.back().offset = newpos.first - params.graph.NodeStart(newNodeIndex);
			result.back().readpos = newpos.second;
			result.back().graphChar = params.graph.NodeSequences(newpos.first);
			result.back().readChar = sequence[newpos.second];
			if (newpos.second == oldpos.second)
			{
				result.back().type = AlignmentResult::TraceMatchType::DELETION;
			}
			else if (newpos.first == oldpos.first && !diagonal)
			{
				result.back().type = AlignmentResult::TraceMatchType::INSERTION;
			}
			else
			{
				assert(diagonal);
				if (Common::characterMatch(sequence[newpos.second], params.graph.NodeSequences(newpos.first)))
				{
					result.back().type = AlignmentResult::TraceMatchType::MATCH;
				}
				else
				{
					result.back().type = AlignmentResult::TraceMatchType::MISMATCH;
				}
			}
		}
		return result;
	}

#ifndef NDEBUG
	void verifyTrace(const std::vector<MatrixPosition>& trace, const std::string& sequence, volatile ScoreType score) const
	{
		volatile ScoreType realscore = 0;
		assert(trace[0].second == 0);
		realscore += Common::characterMatch(sequence[0], params.graph.NodeSequences(trace[0].first)) ? 0 : 1;
		for (size_t i = 1; i < trace.size(); i++)
		{
			auto newpos = trace[i];
			auto oldpos = trace[i-1];
			assert(newpos.second == oldpos.second || newpos.second == oldpos.second+1);
			assert(newpos.second != oldpos.second || newpos.first != oldpos.first);
			auto oldNodeIndex = params.graph.IndexToNode(oldpos.first);
			if (oldpos.first == params.graph.NodeEnd(oldNodeIndex)-1)
			{
				auto newNodeIndex = params.graph.IndexToNode(newpos.first);
				assert(newpos.first == oldpos.first || newpos.first == params.graph.NodeStart(newNodeIndex));
			}
			else
			{
				assert(newpos.first == oldpos.first || newpos.first == oldpos.first+1);
			}
			bool diagonal = true;
			if (newpos.second == oldpos.second) diagonal = false;
			if (newpos.first == oldpos.first)
			{
				auto newNodeIndex = params.graph.IndexToNode(newpos.first);
				if (newpos.second == oldpos.second+1 && params.graph.NodeEnd(newNodeIndex) == params.graph.NodeStart(newNodeIndex)+1 && std::find(params.graph.outNeighbors[newNodeIndex].begin(), params.graph.outNeighbors[newNodeIndex].end(), newNodeIndex) != params.graph.outNeighbors[newNodeIndex].end())
				{
					//one node self-loop, diagonal is valid
				}
				else
				{
					diagonal = false;
				}
			}
			if (!diagonal || !Common::characterMatch(sequence[newpos.second], params.graph.NodeSequences(newpos.first)))
			{
				realscore++;
			}
		}
		// assert(score == realscore);
	}
#endif

};

#endif