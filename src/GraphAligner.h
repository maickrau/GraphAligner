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
#include "GraphAlignerGAFAlignment.h"
#include "GraphAlignerBitvectorBanded.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAligner
{
private:
	using VGAlignment = GraphAlignerVGAlignment<LengthType, ScoreType, Word>;
	using GAFAlignment = GraphAlignerGAFAlignment<LengthType, ScoreType, Word>;
	using BitvectorAligner = GraphAlignerBitvectorBanded<LengthType, ScoreType, Word>;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	using TraceItem = typename Common::TraceItem;
	const Params& params;
	BitvectorAligner bvAligner;
	mutable BufferedWriter logger;
public:

	GraphAligner(const Params& params) :
	params(params),
	bvAligner(params),
	logger()
	{
		if (!params.quietMode) logger = { std::cerr };
	}
	
	AlignmentResult AlignOneWay(const std::string& seq_id, const std::string& sequence, AlignerGraphsizedState& reusableState, size_t DPRestartStride) const
	{
		AlignmentResult result;
		result.readName = seq_id;
		std::string bwSequence = CommonUtils::ReverseComplement(sequence);
		auto fw = fullstartOneWay(seq_id, reusableState, sequence, bwSequence, 0);
		if (!fw.alignmentFailed()) result.alignments.emplace_back(std::move(fw));
		if (DPRestartStride > 0)
		{
			size_t start = 0;
			size_t lastEnd = 0;
			if (result.alignments.size() > 0) lastEnd = result.alignments.back().alignmentEnd;
			while (start < sequence.size())
			{
				start = lastEnd + DPRestartStride;
				if (start >= sequence.size()-1) break;
				auto aln = fullstartOneWay(seq_id, reusableState, sequence, bwSequence, start);
				if (!aln.alignmentFailed())
				{
					assert(aln.alignmentEnd > lastEnd);
					lastEnd = aln.alignmentEnd;
					result.alignments.emplace_back(std::move(aln));
				}
				else
				{
					lastEnd = start;
				}
			}
		}
		return result;
	}

	AlignmentResult AlignClusters(const std::string& seq_id, const std::string& sequence, const std::vector<SeedCluster>& seedClusters, AlignerGraphsizedState& reusableState) const
	{
		assert(params.graph.Finalized());
		AlignmentResult result;
		result.readName = seq_id;
		assert(seedClusters.size() > 0);
		std::string revSequence = CommonUtils::ReverseComplement(sequence);
		std::vector<ScoreType> sliceMaxScores;
		sliceMaxScores.resize(sequence.size() / WordConfiguration<Word>::WordSize + 2, 0);
		for (size_t i = 0; i < seedClusters.size(); i++)
		{
			if (!logger.inputDiscarded()) logger << seq_id << " cluster " << i << "/" << seedClusters.size() << " " << seedClusters[i].clusterGoodness;
			logger << BufferedWriter::Flush;
			result.seedsExtended += 1;
			auto alns = getAlignmentsFromMultiseeds(sequence, revSequence, seedClusters[i].hits, reusableState, sliceMaxScores);
			if (alns.size() == 0) continue;
			for (auto& item : alns)
			{
				assert(!item.alignmentFailed());
				item.seedGoodness = seedClusters[i].clusterGoodness;
				result.alignments.emplace_back(std::move(item));
			}
		}
		assertSetNoRead(seq_id);

		return result;
	}

	void AddAlignment(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment) const
	{
		assert(alignment.trace->trace.size() > 0);
		auto vgAln = VGAlignment::traceToAlignment(seq_id, sequence, alignment.trace->score, alignment.trace->trace, alignment.mappingQuality, 0, false);
		alignment.alignment = vgAln;
		alignment.alignment->set_sequence(sequence.substr(alignment.alignmentStart, alignment.alignmentEnd - alignment.alignmentStart));
		alignment.alignment->set_query_position(alignment.alignmentStart);
	}

	void AddGAFLine(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment, bool cigarMatchMismatchMerge, const bool includeCigar) const
	{
		assert(alignment.trace->trace.size() > 0);
		alignment.GAFline = GAFAlignment::traceToAlignment(seq_id, sequence, *alignment.trace, alignment.alignmentXScore, alignment.mappingQuality, params, cigarMatchMismatchMerge, includeCigar);
	}

	void AddCorrected(AlignmentResult::AlignmentItem& alignment) const
	{
		assert(alignment.trace != nullptr);
		assert(alignment.trace->trace.size() > 0);
		alignment.corrected.reserve(alignment.trace->trace.back().DPposition.seqPos - alignment.trace->trace[0].DPposition.seqPos);
		alignment.corrected = alignment.trace->trace[0].graphCharacter;
		for (size_t i = 1; i < alignment.trace->trace.size(); i++)
		{
			if (!alignment.trace->trace[i-1].nodeSwitch && alignment.trace->trace[i].DPposition.nodeOffset == alignment.trace->trace[i-1].DPposition.nodeOffset && alignment.trace->trace[i].DPposition.node == alignment.trace->trace[i-1].DPposition.node) continue;
			alignment.corrected += alignment.trace->trace[i].graphCharacter;
		}
	}

	std::vector<SeedCluster> clusterSeeds(const std::vector<SeedHit>& seedHits, const size_t seedClusterMinSize) const
	{
		phmap::flat_hash_map<size_t, std::vector<std::tuple<size_t, size_t, size_t>>> seedPoses;
		for (size_t i = 0; i < seedHits.size(); i++)
		{
			int forwardNodeId;
			if (seedHits[i].reverse)
			{
				forwardNodeId = seedHits[i].nodeID * 2 + 1;
			}
			else
			{
				forwardNodeId = seedHits[i].nodeID * 2;
			}
			size_t nodeIndex = params.graph.GetDigraphNode(forwardNodeId, seedHits[i].nodeOffset);
			assert(forwardNodeId == params.graph.BigraphNodeID(nodeIndex));
			assert(seedHits[i].nodeOffset >= params.graph.NodeOffset(nodeIndex));
			assert(seedHits[i].nodeOffset < params.graph.NodeOffset(nodeIndex) + params.graph.NodeLength(nodeIndex));
			assert(params.graph.ChainApproxPos(forwardNodeId) + seedHits[i].nodeOffset > seedHits[i].seqPos);
			seedPoses[params.graph.ChainNumber(forwardNodeId)].emplace_back(i, params.graph.ChainApproxPos(forwardNodeId) + seedHits[i].nodeOffset - seedHits[i].seqPos, nodeIndex);
		}
		std::vector<SeedCluster> clusters;
		for (auto& pair : seedPoses)
		{
			std::sort(pair.second.begin(), pair.second.end(), [](std::tuple<size_t, size_t, size_t> left, std::tuple<size_t, size_t, size_t> right) { return std::get<1>(left) < std::get<1>(right); });
			size_t clusterStart = 0;
			for (size_t i = 1; i <= pair.second.size(); i++)
			{
				assert(i == pair.second.size() || std::get<1>(pair.second[i]) >= std::get<1>(pair.second[i-1]));
				if (i < pair.second.size() && std::get<1>(pair.second[i]) <= std::get<1>(pair.second[i-1]) + 100) continue;
				assert(i > clusterStart);
				if (i - clusterStart < seedClusterMinSize)
				{
					clusterStart = i;
					continue;
				}
				clusters.emplace_back();
				std::sort(pair.second.begin()+clusterStart, pair.second.begin()+i, [&seedHits](std::tuple<size_t, size_t, size_t> left, std::tuple<size_t, size_t, size_t> right) { return seedHits[std::get<0>(left)].seqPos < seedHits[std::get<0>(right)].seqPos; });
				size_t matchingBps = 0;
				int lastEnd = std::numeric_limits<int>::min();
				for (size_t j = clusterStart; j < i; j++)
				{
					size_t seedHitId = std::get<0>(pair.second[j]);
					int thisStart = (int)seedHits[seedHitId].seqPos - (int)seedHits[seedHitId].matchLen + 1;
					int thisEnd = (int)seedHits[seedHitId].seqPos;
					assert(thisEnd >= lastEnd);
					assert(thisEnd > thisStart);
					matchingBps += (thisEnd - std::max(thisStart, lastEnd));
					lastEnd = thisEnd;
				}
				clusters.back().clusterGoodness = matchingBps;
				std::unordered_set<size_t> nodeIdsInSlice;
				size_t lastSlice = 0;
				for (size_t j = clusterStart; j < i; j++)
				{
					size_t seedHitId = std::get<0>(pair.second[j]);
					size_t thisSlice = seedHits[seedHitId].seqPos / WordConfiguration<Word>::WordSize;
					if (thisSlice != lastSlice)
					{
						for (auto id : nodeIdsInSlice)
						{
							clusters.back().hits.emplace_back(lastSlice * WordConfiguration<Word>::WordSize, id);
						}
						lastSlice = thisSlice;
						nodeIdsInSlice.clear();
					}
					nodeIdsInSlice.insert(std::get<2>(pair.second[j]));
				}
				for (auto id : nodeIdsInSlice)
				{
					clusters.back().hits.emplace_back(lastSlice * WordConfiguration<Word>::WordSize, id);
				}
				clusterStart = i;
			}
		}
		std::sort(clusters.begin(), clusters.end(), [](const SeedCluster& left, const SeedCluster& right) { return left.clusterGoodness > right.clusterGoodness; });
		return clusters;
	}

private:

	void fixOverlapTrace(OnewayTrace& trace) const
	{
		fixOverlapTraceStart(trace);
		fixOverlapTraceEnd(trace);
	}

	void fixOverlapTraceEnd(OnewayTrace& trace) const
	{
		size_t fixyStart = trace.trace.size()-1;
		size_t seqAfterHere = 0;
		for (size_t i = trace.trace.size()-1; i > 0; i--)
		{
			if (trace.trace[i].DPposition.nodeOffset != trace.trace[i-1].DPposition.nodeOffset || trace.trace[i].DPposition.node != trace.trace[i-1].DPposition.node || trace.trace[i].nodeSwitch)
			{
				seqAfterHere += 1;
			}
			if (trace.trace[i].DPposition.node == trace.trace[i-1].DPposition.node && !trace.trace[i].nodeSwitch) continue;
			if (seqAfterHere > params.graph.BigraphNodeSize(trace.trace[i-1].DPposition.node) - trace.trace[i-1].DPposition.nodeOffset-1) break;
			fixyStart = i-1;
		}
		if (fixyStart == trace.trace.size()-1) return;
		assert(fixyStart < trace.trace.size());
		std::vector<AlignmentGraph::MatrixPosition> fixyPart;
		fixyPart.reserve(trace.trace.size());
		AlignmentGraph::MatrixPosition graphPos = trace.trace[fixyStart].DPposition;
		for (size_t i = fixyStart+1; i < trace.trace.size(); i++)
		{
			if (trace.trace[i].DPposition.seqPos != trace.trace[i-1].DPposition.seqPos)
			{
				assert(trace.trace[i].DPposition.seqPos == trace.trace[i-1].DPposition.seqPos+1);
				graphPos.seqPos += 1;
			}
			if (trace.trace[i].DPposition.nodeOffset != trace.trace[i-1].DPposition.nodeOffset || trace.trace[i].DPposition.node != trace.trace[i-1].DPposition.node || trace.trace[i].nodeSwitch)
			{
				assert(graphPos.nodeOffset < params.graph.BigraphNodeSize(graphPos.node)-1);
				graphPos.nodeOffset += 1;
			}
			size_t unitigNode = params.graph.GetDigraphNode(graphPos.node, graphPos.nodeOffset);
			size_t nodeOffset = params.graph.NodeOffset(unitigNode);
			if (params.graph.NodeSequences(unitigNode, graphPos.nodeOffset - nodeOffset) != trace.trace[i-1].graphCharacter)
			{
				// imperfect overlap, can't be correctly fixed
				fixyPart.clear();
				break;
			}
			fixyPart.push_back(graphPos);
		}
		if (fixyPart.size() != trace.trace.size()-1-fixyStart) return;
		for (size_t i = 0; i < fixyPart.size(); i++)
		{
			trace.trace[fixyStart+1+i].DPposition = fixyPart[i];
		}
	}

	void fixOverlapTraceStart(OnewayTrace& trace) const
	{
		size_t fixyStart = 0;
		size_t seqBeforeHere = 0;
		for (size_t i = 1; i < trace.trace.size(); i++)
		{
			if (trace.trace[i].DPposition.nodeOffset != trace.trace[i-1].DPposition.nodeOffset || trace.trace[i].DPposition.node != trace.trace[i-1].DPposition.node || trace.trace[i].nodeSwitch)
			{
				seqBeforeHere += 1;
			}
			if (trace.trace[i].DPposition.nodeOffset < seqBeforeHere) break;
			if (trace.trace[i].DPposition.node == trace.trace[i-1].DPposition.node && !trace.trace[i].nodeSwitch) continue;
			fixyStart = i;
		}
		if (fixyStart == 0) return;
		assert(fixyStart < trace.trace.size());
		std::vector<AlignmentGraph::MatrixPosition> fixyPart;
		fixyPart.reserve(trace.trace.size());
		AlignmentGraph::MatrixPosition graphPos = trace.trace[fixyStart].DPposition;
		for (size_t i = fixyStart; i > 0; i--)
		{
			if (trace.trace[i].DPposition.seqPos != trace.trace[i-1].DPposition.seqPos)
			{
				assert(trace.trace[i].DPposition.seqPos == trace.trace[i-1].DPposition.seqPos+1);
				graphPos.seqPos -= 1;
			}
			if (trace.trace[i].DPposition.nodeOffset != trace.trace[i-1].DPposition.nodeOffset || trace.trace[i].nodeSwitch || trace.trace[i].DPposition.node != trace.trace[i-1].DPposition.node)
			{
				assert(graphPos.nodeOffset != 0);
				graphPos.nodeOffset -= 1;
			}
			size_t unitigNode = params.graph.GetDigraphNode(graphPos.node, graphPos.nodeOffset);
			size_t nodeOffset = params.graph.NodeOffset(unitigNode);
			if (params.graph.NodeSequences(unitigNode, graphPos.nodeOffset - nodeOffset) != trace.trace[i-1].graphCharacter)
			{
				// imperfect overlap, can't be correctly fixed
				fixyPart.clear();
				break;
			}
			fixyPart.push_back(graphPos);
		}
		if (fixyPart.size() != fixyStart) return;
		for (size_t i = 0; i < fixyPart.size(); i++)
		{
			trace.trace[fixyStart-1-i].DPposition = fixyPart[i];
		}
	}

	OnewayTrace mergeTraces(OnewayTrace&& bwTrace, OnewayTrace&& fwTrace) const
	{
		assert(!fwTrace.failed() || !bwTrace.failed());
		OnewayTrace mergedTrace;
		if (!bwTrace.failed())
		{
			assert(bwTrace.trace.size() > 0);
			mergedTrace = std::move(bwTrace);
			if (!fwTrace.failed())
			{
				assert(mergedTrace.trace.size() > 0);
				assert(fwTrace.trace.size() > 0);
				assert(mergedTrace.trace.back().DPposition == fwTrace.trace[0].DPposition);
				if (!Common::characterMatch(fwTrace.trace[0].sequenceCharacter, fwTrace.trace[0].graphCharacter))
				{
					assert(bwTrace.score >= 1);
					assert(fwTrace.score >= 1);
					fwTrace.score -= 1;
				}
				mergedTrace.trace.pop_back();
				mergedTrace.trace.insert(mergedTrace.trace.end(), fwTrace.trace.begin(), fwTrace.trace.end());
				mergedTrace.score += fwTrace.score;
			}
		}
		else if (!fwTrace.failed())
		{
			mergedTrace = std::move(fwTrace);
		}
		return std::move(mergedTrace);
	}

	OnewayTrace clipAndAddBackwardTrace(OnewayTrace&& fwTrace, AlignerGraphsizedState& reusableState, const std::string& fwSequence, const std::string& bwSequence, size_t offset) const
	{
		clipTraceStart(fwTrace);
		if (fwTrace.trace.size() == 0) return std::move(fwTrace);
		fixForwardTraceSeqPos(fwTrace, offset, fwSequence);
		OnewayTrace bwTrace = OnewayTrace::TraceFailed();
		if (fwTrace.trace[0].DPposition.seqPos != 0)
		{
			std::string_view backwardPart { bwSequence.data() + bwSequence.size() - fwTrace.trace[0].DPposition.seqPos, fwTrace.trace[0].DPposition.seqPos };
			auto reversePos = params.graph.GetReversePosition(fwTrace.trace[0].DPposition.node, fwTrace.trace[0].DPposition.nodeOffset);
			bwTrace = bvAligner.getReverseTraceFromSeed(backwardPart, reversePos.first, reversePos.second, params.Xdropcutoff, reusableState);
			if (!bwTrace.failed())
			{
#ifdef EXTRACORRECTNESSASSERTIONS
				std::reverse(bwTrace.trace.begin(), bwTrace.trace.end());
				verifyTrace(bwTrace.trace, backwardPart, bwTrace.score);
				std::reverse(bwTrace.trace.begin(), bwTrace.trace.end());
#endif
				fixReverseTraceSeqPosAndOrder(bwTrace, fwTrace.trace[0].DPposition.seqPos-1, fwSequence, bwSequence);
			}
		}

		OnewayTrace mergedTrace = mergeTraces(std::move(bwTrace), std::move(fwTrace));
#ifndef NDEBUG
		if (mergedTrace.trace.size() > 0) verifyFixedTrace(mergedTrace.trace, fwSequence, mergedTrace.score);
#endif
		return std::move(mergedTrace);
	}

	AlignmentResult::AlignmentItem fullstartOneWay(const std::string& seq_id, AlignerGraphsizedState& reusableState, const std::string& fwSequence, const std::string& bwSequence, size_t offset) const
	{
		auto timeStart = std::chrono::system_clock::now();
		assert(params.graph.Finalized());
		std::string_view fwView { fwSequence.data() + offset, fwSequence.size() - offset };
		auto fwTrace = getBacktraceFullStart(fwView, reusableState);
		auto timeEnd = std::chrono::system_clock::now();
		//failed alignment, don't output
		if (fwTrace.score == std::numeric_limits<ScoreType>::max()) return AlignmentResult::AlignmentItem {};
		if (fwTrace.trace.size() == 0) return AlignmentResult::AlignmentItem {};
#ifndef NDEBUG
		if (fwTrace.trace.size() > 0) verifyTrace(fwTrace.trace, fwView, fwTrace.score);
#endif
		OnewayTrace mergedTrace = clipAndAddBackwardTrace(std::move(fwTrace), reusableState, fwSequence, bwSequence, offset);
		if (mergedTrace.trace.size() == 0) return AlignmentResult::AlignmentItem {};
		fixOverlapTrace(mergedTrace);

		AlignmentResult::AlignmentItem alnItem { std::move(mergedTrace), 0, std::numeric_limits<size_t>::max() };

		alnItem.alignmentScore = alnItem.trace->score;
		alnItem.alignmentStart = alnItem.trace->trace[0].DPposition.seqPos;
		alnItem.alignmentEnd = alnItem.trace->trace.back().DPposition.seqPos + 1;
		timeEnd = std::chrono::system_clock::now();
		auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
		alnItem.elapsedMilliseconds = time;
		alnItem.alignmentXScore = (ScoreType)alnItem.alignmentLength()*100 - params.XscoreErrorCost * (ScoreType)alnItem.alignmentScore;
		return alnItem;
	}

	// bool exactAlignmentPart(const AlignmentResult::AlignmentItem& aln, const SeedHit& seedHit) const
	// {
	// 	assert(aln.trace != nullptr);
	// 	const std::vector<TraceItem>& trace = aln.trace->trace;
	// 	assert(trace.size() > 0);
	// 	assert(trace.back().DPposition.seqPos > trace[0].DPposition.seqPos);
	// 	if (trace.back().DPposition.seqPos < seedHit.seqPos) return false;
	// 	if (trace[0].DPposition.seqPos > seedHit.seqPos) return false;
	// 	size_t high = trace.size();
	// 	size_t low = 0;
	// 	size_t mid = (seedHit.seqPos - trace[0].DPposition.seqPos) / (trace.back().DPposition.seqPos - trace[0].DPposition.seqPos);
	// 	while (trace[mid].DPposition.seqPos != seedHit.seqPos)
	// 	{
	// 		if (trace[mid].DPposition.seqPos < seedHit.seqPos)
	// 		{
	// 			low = mid;
	// 			mid = (high + low) / 2;
	// 			if (mid == low) mid += 1;
	// 			assert(mid < trace.size());
	// 		}
	// 		if (trace[mid].DPposition.seqPos > seedHit.seqPos)
	// 		{
	// 			high = mid;
	// 			mid = (high + low) / 2;
	// 			assert(mid < trace.size());
	// 		}
	// 		assert(low < mid);
	// 		assert(mid < high);
	// 	}
	// 	assert(mid < trace.size());
	// 	assert(trace[mid].DPposition.seqPos == seedHit.seqPos);
	// 	size_t down = mid;
	// 	size_t compareNode = seedHit.nodeID * 2;
	// 	if (seedHit.reverse) compareNode += 1;
	// 	while (trace[down].DPposition.seqPos == seedHit.seqPos)
	// 	{
	// 		if (compareNode == trace[down].DPposition.node && seedHit.nodeOffset == trace[down].DPposition.nodeOffset)
	// 		{
	// 			return true;
	// 		}
	// 		if (down == 0) break;
	// 		down -= 1;
	// 	}
	// 	size_t up = mid;
	// 	while (trace[up].DPposition.seqPos == seedHit.seqPos)
	// 	{
	// 		if (compareNode == trace[up].DPposition.node && seedHit.nodeOffset == trace[up].DPposition.nodeOffset)
	// 		{
	// 			return true;
	// 		}
	// 		up += 1;
	// 		if (up == trace.size()) break;
	// 	}
	// 	return false;
	// }

	OnewayTrace getBacktraceFullStart(const std::string& sequence, AlignerGraphsizedState& reusableState) const
	{
		std::string_view seq { sequence.data(), sequence.size() };
		return getBacktraceFullStart(seq, reusableState);
	}

	OnewayTrace getBacktraceFullStart(const std::string_view& seq, AlignerGraphsizedState& reusableState) const
	{
		return bvAligner.getBacktraceFullStart(seq, params.Xdropcutoff, reusableState);
	}

	std::vector<OnewayTrace> getMultiseedTraces(const std::string& sequence, const std::string& revSequence, const std::vector<ProcessedSeedHit>& seedHits, AlignerGraphsizedState& reusableState, std::vector<ScoreType>& sliceMaxScores) const
	{
		return bvAligner.getMultiseedTraces(sequence, revSequence, seedHits, reusableState, sliceMaxScores);
	}

	// Trace getTwoDirectionalTrace(const std::string& sequence, const std::string& revSequence, SeedHit seedHit, AlignerGraphsizedState& reusableState) const
	// {
	// 	assert(seedHit.seqPos >= 0);
	// 	assert(seedHit.seqPos < sequence.size());
	// 	int forwardNodeId;
	// 	int backwardNodeId;
	// 	if (seedHit.reverse)
	// 	{
	// 		forwardNodeId = seedHit.nodeID * 2 + 1;
	// 		backwardNodeId = seedHit.nodeID * 2;
	// 	}
	// 	else
	// 	{
	// 		forwardNodeId = seedHit.nodeID * 2;
	// 		backwardNodeId = seedHit.nodeID * 2 + 1;
	// 	}
	// 	Trace result;
	// 	result.backward.score = std::numeric_limits<ScoreType>::max();
	// 	result.forward.score = std::numeric_limits<ScoreType>::max();
	// 	if (seedHit.seqPos > 0)
	// 	{
	// 		std::string_view backwardPart { revSequence.data() + revSequence.size() - seedHit.seqPos, seedHit.seqPos };
	// 		auto reversePos = params.graph.GetReversePosition(forwardNodeId, seedHit.nodeOffset);
	// 		assert(reversePos.first == backwardNodeId);
	// 		result.backward = bvAligner.getReverseTraceFromSeed(backwardPart, backwardNodeId, reversePos.second, params.Xdropcutoff, reusableState);
	// 	}
	// 	if (seedHit.seqPos < sequence.size()-1)
	// 	{
	// 		std::string_view forwardPart { sequence.data() + seedHit.seqPos + 1, sequence.size() - seedHit.seqPos - 1 };
	// 		size_t offset = seedHit.nodeOffset;
	// 		result.forward = bvAligner.getReverseTraceFromSeed(forwardPart, forwardNodeId, offset, params.Xdropcutoff, reusableState);
	// 	}

	// 	if (!result.backward.failed())
	// 	{
	// 		auto reversePos = params.graph.GetReversePosition(forwardNodeId, seedHit.nodeOffset);
	// 		assert(result.backward.trace.back().DPposition.seqPos == (size_t)-1 && params.graph.BigraphNodeID(result.backward.trace.back().DPposition.node) == backwardNodeId && params.graph.NodeOffset(result.backward.trace.back().DPposition.node) + result.backward.trace.back().DPposition.nodeOffset == reversePos.second);
	// 		std::reverse(result.backward.trace.begin(), result.backward.trace.end());
	// 	}
	// 	if (!result.forward.failed())
	// 	{
	// 		assert(result.forward.trace.back().DPposition.seqPos == (size_t)-1 && params.graph.BigraphNodeID(result.forward.trace.back().DPposition.node) == forwardNodeId && params.graph.NodeOffset(result.forward.trace.back().DPposition.node) + result.forward.trace.back().DPposition.nodeOffset == seedHit.nodeOffset);
	// 		std::reverse(result.forward.trace.begin(), result.forward.trace.end());
	// 	}
	// 	return result;
	// }

	void fixForwardTraceSeqPos(OnewayTrace& trace, LengthType start, const std::string& sequence) const
	{
		if (trace.trace.size() == 0) return;
		trace.score = 0;
		for (size_t i = 0; i < trace.trace.size(); i++)
		{
			trace.trace[i].DPposition.seqPos += start;
			auto nodeIndex = trace.trace[i].DPposition.node;
			trace.trace[i].DPposition.node = params.graph.BigraphNodeID(nodeIndex);
			trace.trace[i].DPposition.nodeOffset += params.graph.NodeOffset(nodeIndex);
			assert(trace.trace[i].DPposition.seqPos >= 0);
			assert(trace.trace[i].DPposition.seqPos < sequence.size());
			if (trace.trace[i].DPposition.seqPos == start-1)
			{
				assert(trace.trace[i].sequenceCharacter == '-');
				trace.trace[i].sequenceCharacter = sequence[trace.trace[i].DPposition.seqPos];
			}
			assert(trace.trace[i].sequenceCharacter == sequence[trace.trace[i].DPposition.seqPos]);
			if (i == 0 && !Common::characterMatch(trace.trace[0].sequenceCharacter, trace.trace[0].graphCharacter)) trace.score += 1;
			if (i > 0)
			{
				if (trace.trace[i].DPposition.seqPos == trace.trace[i-1].DPposition.seqPos)
				{
					trace.score += 1;
				}
				else if (!trace.trace[i-1].nodeSwitch && trace.trace[i].DPposition.node == trace.trace[i-1].DPposition.node && trace.trace[i].DPposition.nodeOffset == trace.trace[i-1].DPposition.nodeOffset)
				{
					trace.score += 1;
				}
				else if (!Common::characterMatch(trace.trace[i].sequenceCharacter, trace.trace[i].graphCharacter))
				{
					trace.score += 1;
				}
			}
			assert(trace.trace[i].DPposition.seqPos < sequence.size());
			assert(i == 0 || trace.trace[i].DPposition.seqPos == trace.trace[0].DPposition.seqPos || trace.trace[i].sequenceCharacter == sequence[trace.trace[i].DPposition.seqPos]);
		}
	}

	//end is the (fw) index of the last alignable base pair, not one beyond
	void fixReverseTraceSeqPosAndOrder(OnewayTrace& trace, LengthType end, const std::string& sequence, const std::string& revSequence) const
	{
		if (trace.trace.size() == 0) return;
		std::vector<bool> diagonalIndices;
		diagonalIndices.resize(trace.trace.size(), false);
		for (size_t i = 0; i < trace.trace.size(); i++)
		{
			bool diagonal = false;
			if (i == trace.trace.size()-1)
			{
				diagonal = true;
			}
			else
			{
				if (trace.trace[i].DPposition.seqPos != trace.trace[i+1].DPposition.seqPos)
				{
					assert(trace.trace[i].DPposition.seqPos == trace.trace[i+1].DPposition.seqPos+1);
					if (trace.trace[i].DPposition.node != trace.trace[i+1].DPposition.node || trace.trace[i].DPposition.nodeOffset != trace.trace[i+1].DPposition.nodeOffset || trace.trace[i+1].nodeSwitch)
					{
						diagonal = true;
					}
				}
			}
			diagonalIndices[i] = diagonal;
		}
		assert(diagonalIndices[0]);
		assert(diagonalIndices.back());
		for (size_t i = 0; i < trace.trace.size() - 1; i++)
		{
			trace.trace[i].nodeSwitch = trace.trace[i+1].nodeSwitch;
		}
		trace.trace.back().nodeSwitch = false;
		for (size_t i = 0; i < trace.trace.size(); i++)
		{
			assert(trace.trace[i].DPposition.seqPos <= end || trace.trace[i].DPposition.seqPos == (size_t)-1);
			trace.trace[i].DPposition.seqPos = end - trace.trace[i].DPposition.seqPos;
			size_t offset = params.graph.NodeOffset(trace.trace[i].DPposition.node) + trace.trace[i].DPposition.nodeOffset;
			auto reversePos = params.graph.GetReversePosition(params.graph.BigraphNodeID(trace.trace[i].DPposition.node), offset);
			assert(reversePos.second < params.graph.BigraphNodeSize(params.graph.BigraphNodeID(trace.trace[i].DPposition.node)));
			trace.trace[i].DPposition.node = reversePos.first;
			trace.trace[i].DPposition.nodeOffset = reversePos.second;
			assert(trace.trace[i].DPposition.seqPos < sequence.size());
			trace.trace[i].sequenceCharacter = sequence[trace.trace[i].DPposition.seqPos];
			trace.trace[i].graphCharacter = CommonUtils::Complement(trace.trace[i].graphCharacter);
		}
		trace.score = 0;
		size_t lastDiagonalPos;
		bool lastNodeSwitch;
		for (size_t i = 0; i < trace.trace.size(); i++)
		{
			if (diagonalIndices[i])
			{
				lastDiagonalPos = i;
				lastNodeSwitch = trace.trace[i].nodeSwitch;
				if (!Common::characterMatch(trace.trace[i].sequenceCharacter, trace.trace[i].graphCharacter)) trace.score += 1;
			}
			else
			{
				assert(i < trace.trace.size()-1);
				assert((trace.trace[i].DPposition.seqPos != trace.trace[i+1].DPposition.seqPos && trace.trace[i].DPposition.node == trace.trace[i+1].DPposition.node && trace.trace[i].DPposition.nodeOffset == trace.trace[i+1].DPposition.nodeOffset && !trace.trace[i].nodeSwitch) || (trace.trace[i].DPposition.seqPos == trace.trace[i+1].DPposition.seqPos && (trace.trace[i].DPposition.node != trace.trace[i+1].DPposition.node || trace.trace[i].DPposition.nodeOffset != trace.trace[i+1].DPposition.nodeOffset || trace.trace[i].nodeSwitch)));
				if (trace.trace[i].DPposition.seqPos != trace.trace[i+1].DPposition.seqPos)
				{
					trace.trace[i].DPposition.node = trace.trace[lastDiagonalPos].DPposition.node;
					trace.trace[i].DPposition.nodeOffset = trace.trace[lastDiagonalPos].DPposition.nodeOffset;
					trace.trace[lastDiagonalPos].nodeSwitch = false;
					if (!diagonalIndices[i+1])
					{
						trace.trace[i].nodeSwitch = false;
					}
					else
					{
						trace.trace[i].nodeSwitch = lastNodeSwitch;
					}
				}
				else
				{
					trace.trace[i].DPposition.seqPos = trace.trace[lastDiagonalPos].DPposition.seqPos;
				}
				trace.score += 1;
			}
		}
		for (size_t i = 1; i < trace.trace.size(); i++)
		{
			assert(trace.trace[i].DPposition.seqPos == trace.trace[i-1].DPposition.seqPos || trace.trace[i].DPposition.seqPos == trace.trace[i-1].DPposition.seqPos+1);
		}
	}

	void fixTraceConsecutiveIndels(OnewayTrace& trace) const
	{
		size_t skipped = 0;
		for (size_t i = 1; i < trace.trace.size()-1; i++)
		{
			size_t previousI = i - skipped - 1;
			bool prevDel = false;
			bool prevIns = false;
			bool nextDel = false;
			bool nextIns = false;
			if (trace.trace[i].DPposition.seqPos == trace.trace[previousI].DPposition.seqPos)
			{
				prevDel = true;
				assert(trace.trace[previousI].nodeSwitch || trace.trace[i].DPposition.node != trace.trace[previousI].DPposition.node || trace.trace[i].DPposition.nodeOffset != trace.trace[previousI].DPposition.nodeOffset);
			}
			if (!trace.trace[previousI].nodeSwitch && trace.trace[previousI].DPposition.node == trace.trace[i].DPposition.node && trace.trace[previousI].DPposition.nodeOffset == trace.trace[i].DPposition.nodeOffset)
			{
				prevIns = true;
				assert(trace.trace[previousI].DPposition.seqPos != trace.trace[i].DPposition.seqPos);
			}
			assert(!prevDel || !prevIns);
			if (trace.trace[i].DPposition.seqPos == trace.trace[i+1].DPposition.seqPos)
			{
				nextDel = true;
				assert(trace.trace[i].nodeSwitch || trace.trace[i].DPposition.node != trace.trace[i+1].DPposition.node || trace.trace[i].DPposition.nodeOffset != trace.trace[i+1].DPposition.nodeOffset);
			}
			if (!trace.trace[i].nodeSwitch && trace.trace[i].DPposition.node == trace.trace[i+1].DPposition.node && trace.trace[i].DPposition.nodeOffset == trace.trace[i+1].DPposition.nodeOffset)
			{
				nextIns = true;
				assert(trace.trace[i].DPposition.seqPos != trace.trace[i+1].DPposition.seqPos);
			}
			assert(!nextIns || !nextDel);
			if (prevDel && nextIns)
			{
				skipped += 1;
				if (trace.trace[i].nodeSwitch) trace.trace[previousI].nodeSwitch = true;
			}
			else if (prevIns && nextDel)
			{
				skipped += 1;
				if (trace.trace[i].nodeSwitch) trace.trace[previousI].nodeSwitch = true;
			}
			else
			{
				if (skipped > 0)
				{
					trace.trace[i-skipped] = trace.trace[i];
				}
			}
		}
		if (skipped > 0)
		{
			trace.trace[trace.trace.size()-1-skipped] = trace.trace[trace.trace.size()-1];
			for (size_t i = 0; i < skipped; i++)
			{
				trace.trace.pop_back();
			}
		}
	}

	std::vector<AlignmentResult::AlignmentItem> getAlignmentsFromMultiseeds(const std::string& sequence, const std::string& revSequence, const std::vector<ProcessedSeedHit>& seedHits, AlignerGraphsizedState& reusableState, std::vector<ScoreType>& sliceMaxScores) const
	{
		auto traces = getMultiseedTraces(sequence, revSequence, seedHits, reusableState, sliceMaxScores);
		std::vector<AlignmentResult::AlignmentItem> result;
		for (size_t i = 0; i < traces.size(); i++)
		{
			assert(!traces[i].failed());
			fixTraceConsecutiveIndels(traces[i]);
			fixReverseTraceSeqPosAndOrder(traces[i], sequence.size()-1, sequence, revSequence);
			fixOverlapTrace(traces[i]);
			ScoreType alignmentXScore = (ScoreType)(traces[i].trace.back().DPposition.seqPos - traces[i].trace[0].DPposition.seqPos + 1)*100 - params.XscoreErrorCost * (ScoreType)traces[i].score;
			if (alignmentXScore <= 0) continue;
			result.emplace_back(std::move(traces[i]), 0, std::numeric_limits<size_t>::max());
			LengthType seqstart = 0;
			LengthType seqend = 0;
			assert(result.back().trace->trace.size() > 0);
			seqstart = result.back().trace->trace[0].DPposition.seqPos;
			seqend = result.back().trace->trace.back().DPposition.seqPos;
			assert(seqend < sequence.size());
			result.back().alignmentScore = result.back().trace->score;
			result.back().alignmentStart = seqstart;
			result.back().alignmentEnd = seqend + 1;
			result.back().alignmentXScore = (ScoreType)result.back().alignmentLength()*100 - params.XscoreErrorCost * (ScoreType)result.back().alignmentScore;
			assert(result.back().alignmentXScore == alignmentXScore);
			result.back().alignmentXScore /= 100.0;
		}
		return result;
	}

// 	AlignmentResult::AlignmentItem getAlignmentFromSeed(const std::string& seq_id, const std::string& sequence, const std::string& revSequence, SeedHit seedHit, AlignerGraphsizedState& reusableState) const
// 	{
// 		assert(params.graph.Finalized());
// 		auto timeStart = std::chrono::system_clock::now();

// 		auto trace = getTwoDirectionalTrace(sequence, revSequence, seedHit, reusableState);

// #ifndef NDEBUG
// 		if (trace.forward.trace.size() > 0) verifyTrace(trace.forward.trace, sequence, trace.forward.score);
// 		if (trace.backward.trace.size() > 0) verifyTrace(trace.backward.trace, sequence, trace.backward.score);
// #endif
// 		fixReverseTraceSeqPosAndOrder(trace.backward, seedHit.seqPos-1, sequence);
// 		fixForwardTraceSeqPos(trace.forward, seedHit.seqPos+1, sequence);

// 		//failed alignment, don't output
// 		if (trace.forward.failed() && trace.backward.failed())
// 		{
// 			return emptyAlignment(0, 0);
// 		}

// 		// auto traceVector = getTraceInfo(sequence, trace.backward.trace, trace.forward.trace);

// 		assert(trace.backward.trace.size() > 0 || trace.backward.failed());
// 		auto mergedTrace = std::move(trace.backward);
// 		if (mergedTrace.failed())
// 		{
// 			mergedTrace = std::move(trace.forward);
// 		}
// 		else if (!mergedTrace.failed() && !trace.forward.failed())
// 		{
// 			assert(mergedTrace.trace.size() > 0);
// 			assert(mergedTrace.trace.back().DPposition == trace.forward.trace[0].DPposition);
// 			mergedTrace.trace.pop_back();
// 			mergedTrace.trace.insert(mergedTrace.trace.end(), trace.forward.trace.begin(), trace.forward.trace.end());
// 			mergedTrace.score += trace.forward.score;
// 		}
// 		else if (!trace.forward.failed())
// 		{
// 			assert(trace.forward.trace.size() > 0);
// 			mergedTrace.trace.insert(mergedTrace.trace.end(), trace.forward.trace.begin(), trace.forward.trace.end());
// 			mergedTrace.score += trace.forward.score;
// 		}

// 		AlignmentResult::AlignmentItem result { std::move(mergedTrace), 0, std::numeric_limits<size_t>::max() };

// 		LengthType seqstart = 0;
// 		LengthType seqend = 0;
// 		assert(result.trace->trace.size() > 0);
// 		seqstart = result.trace->trace[0].DPposition.seqPos;
// 		seqend = result.trace->trace.back().DPposition.seqPos;
// 		assert(seqend < sequence.size());
// 		// result.trace = traceVector;
// 		result.alignmentScore = result.trace->score;
// 		result.alignmentStart = seqstart;
// 		result.alignmentEnd = seqend + 1;
// 		auto timeEnd = std::chrono::system_clock::now();
// 		size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
// 		result.elapsedMilliseconds = time;
// 		result.alignmentXScore = (ScoreType)result.alignmentLength()*100 - params.XscoreErrorCost * (ScoreType)result.alignmentScore;
// 		result.alignmentXScore /= 100.0;
// 		return result;
// 	}

	static AlignmentResult::AlignmentItem emptyAlignment(size_t elapsedMilliseconds, size_t cellsProcessed)
	{
		AlignmentResult::AlignmentItem result;
		result.cellsProcessed = cellsProcessed;
		result.elapsedMilliseconds = elapsedMilliseconds;
		return result;
	}

	void clipTraceStart(OnewayTrace& trace) const
	{
		if (trace.trace.size() == 0) return;
		ScoreType maxX = std::numeric_limits<ScoreType>::min();
		size_t maxXIndex = std::numeric_limits<size_t>::max();
		size_t score = 0;
		size_t maxXScore = 0;
		for (size_t i = trace.trace.size()-1; i < trace.trace.size(); i--)
		{
			bool posSwitch = (i == trace.trace.size()-1) ||  trace.trace[i].nodeSwitch || trace.trace[i].DPposition.nodeOffset != trace.trace[i+1].DPposition.nodeOffset || trace.trace[i].DPposition.node != trace.trace[i+1].DPposition.node;
			if (i != trace.trace.size()-1 && trace.trace[i+1].DPposition.seqPos == trace.trace[i].DPposition.seqPos)
			{
				score += 1;
			}
			else if (i != trace.trace.size()-1 && !posSwitch)
			{
				score += 1;
			}
			else if (!Common::characterMatch(trace.trace[i].sequenceCharacter, trace.trace[i].graphCharacter))
			{
				score += 1;
			}
			ScoreType Xhere = (trace.trace.back().DPposition.seqPos - trace.trace[i].DPposition.seqPos + 1)*100 - score * params.XscoreErrorCost;
			if (Xhere > maxX)
			{
				maxX = Xhere;
				maxXIndex = i;
				maxXScore = score;
			}
		}
		assert(maxXIndex < trace.trace.size());
		if (maxXIndex > 0)
		{
			trace.score = maxXScore;
			trace.trace.erase(trace.trace.begin(), trace.trace.begin() + maxXIndex);
		}
	}

#ifndef NDEBUG
	void verifyTrace(const std::vector<TraceItem>& trace, const std::string& sequence, ScoreType score) const
	{
		std::string_view view { sequence.data(), sequence.size() };
		verifyTrace(trace, view, score);
	}
	void verifyTrace(const std::vector<TraceItem>& trace, const std::string_view& sequence, ScoreType score) const
	{
		for (size_t i = 1; i < trace.size(); i++)
		{
			auto newpos = trace[i].DPposition;
			auto oldpos = trace[i-1].DPposition;
			auto oldNodeIndex = oldpos.node;
			auto newNodeIndex = newpos.node;
			assert(newpos.seqPos != oldpos.seqPos || newpos.node != oldpos.node || newpos.nodeOffset != oldpos.nodeOffset);
			if (oldNodeIndex == newNodeIndex)
			{
				assert(((newpos.nodeOffset == oldpos.nodeOffset || newpos.nodeOffset == oldpos.nodeOffset+1) && (newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1)) || (oldpos.nodeOffset == params.graph.NodeLength(newNodeIndex)-1 && newpos.nodeOffset == 0 && (newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1)));
				continue;
			}
			assert(oldNodeIndex != newNodeIndex);
			assert(std::find(params.graph.OutNeighbors(oldpos.node).begin(), params.graph.OutNeighbors(oldpos.node).end(), newpos.node) != params.graph.OutNeighbors(oldpos.node).end());
			assert(newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1);
			assert(oldpos.nodeOffset == params.graph.NodeLength(oldNodeIndex)-1);
			assert(newpos.nodeOffset == 0);
		}
	}
	void verifyFixedTrace(const std::vector<TraceItem>& trace, const std::string& sequence, ScoreType score) const
	{
		if (trace.size() == 0) return;
		assert(trace[0].DPposition.seqPos < sequence.size());
		assert(trace.back().DPposition.seqPos < sequence.size());
		ScoreType foundScore = 0;
		if (!Common::characterMatch(trace[0].sequenceCharacter, trace[0].graphCharacter)) foundScore += 1;
		for (size_t i = 1; i < trace.size(); i++)
		{
			assert(trace[i].DPposition.seqPos < sequence.size());
			auto newpos = trace[i].DPposition;
			auto oldpos = trace[i-1].DPposition;

			if (newpos.seqPos == oldpos.seqPos)
			{
				foundScore += 1;
			}
			else if (!trace[i-1].nodeSwitch && newpos.node == oldpos.node && newpos.nodeOffset == oldpos.nodeOffset)
			{
				foundScore += 1;
			}
			else if (!Common::characterMatch(trace[i].sequenceCharacter, trace[i].graphCharacter))
			{
				foundScore += 1;
			}

			assert(newpos.nodeOffset < params.graph.BigraphNodeSize(newpos.node));
			size_t nodeIndex = params.graph.GetDigraphNode(newpos.node, newpos.nodeOffset);
			assert(params.graph.NodeOffset(nodeIndex) <= newpos.nodeOffset);
			assert(params.graph.NodeOffset(nodeIndex) + params.graph.NodeLength(nodeIndex) > newpos.nodeOffset);
			size_t offsetInNode = newpos.nodeOffset - params.graph.NodeOffset(nodeIndex);
			assert(offsetInNode < params.graph.NodeLength(nodeIndex));
			newpos.node = nodeIndex;
			newpos.nodeOffset = offsetInNode;

			assert(oldpos.nodeOffset < params.graph.BigraphNodeSize(oldpos.node));
			nodeIndex = params.graph.GetDigraphNode(oldpos.node, oldpos.nodeOffset);
			assert(params.graph.NodeOffset(nodeIndex) <= oldpos.nodeOffset);
			assert(params.graph.NodeOffset(nodeIndex) + params.graph.NodeLength(nodeIndex) > oldpos.nodeOffset);
			offsetInNode = oldpos.nodeOffset - params.graph.NodeOffset(nodeIndex);
			assert(offsetInNode < params.graph.NodeLength(nodeIndex));
			oldpos.node = nodeIndex;
			oldpos.nodeOffset = offsetInNode;

			auto oldNodeIndex = oldpos.node;
			auto newNodeIndex = newpos.node;
			assert(newpos.seqPos != oldpos.seqPos || newpos.node != oldpos.node || newpos.nodeOffset != oldpos.nodeOffset);
			if (oldNodeIndex == newNodeIndex)
			{
				// this check doesn't allow a backwards alignment through a self-looping node with overlap, which is a thing that happens
				// problem: the original backwards trace is fine, aligning to node of size |n| oldpos is at |n|-1 and newpos at 0
				// but when it gets oriented forwards, overlap x, oldpos is at 0 and newpos is at n-x-1 which this check doesn't allow
				// so disable for now
				// assert((newpos.nodeOffset >= oldpos.nodeOffset && newpos.seqPos >= oldpos.seqPos) || (oldpos.nodeOffset == params.graph.NodeLength(newNodeIndex)-1 && newpos.nodeOffset == 0 && (newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1)));
				continue;
			}
			assert(oldNodeIndex != newNodeIndex);
			assert(newpos.seqPos == oldpos.seqPos || newpos.seqPos == oldpos.seqPos+1);
			bool foundSimple = std::find(params.graph.OutNeighbors(oldpos.node).begin(), params.graph.OutNeighbors(oldpos.node).end(), newpos.node) != params.graph.OutNeighbors(oldpos.node).end();
			if (foundSimple)
			{
				assert(newpos.nodeOffset == 0);
				assert(oldpos.nodeOffset == params.graph.NodeLength(oldNodeIndex)-1);
			}
			if (!foundSimple)
			{
				auto revOldNode = (trace[i-1].DPposition.node % 2 == 0) ? (trace[i-1].DPposition.node + 1) : (trace[i-1].DPposition.node - 1);
				auto revNewNode = (trace[i].DPposition.node % 2 == 0) ? (trace[i].DPposition.node + 1) : (trace[i].DPposition.node - 1);
				auto revOldOffset = params.graph.BigraphNodeSize(trace[i-1].DPposition.node) - trace[i-1].DPposition.nodeOffset - 1;
				auto revNewOffset = params.graph.BigraphNodeSize(trace[i].DPposition.node) - trace[i].DPposition.nodeOffset - 1;
				auto revOldNodeIndex = params.graph.GetDigraphNode(revOldNode, revOldOffset);
				auto revNewNodeIndex = params.graph.GetDigraphNode(revNewNode, revNewOffset);
				auto revOldNodeOffset = revOldOffset - params.graph.NodeOffset(revOldNodeIndex);
				auto revNewNodeOffset = revNewOffset - params.graph.NodeOffset(revNewNodeIndex);
				bool foundReverse = std::find(params.graph.OutNeighbors(revNewNodeIndex).begin(), params.graph.OutNeighbors(revNewNodeIndex).end(), revOldNodeIndex) != params.graph.OutNeighbors(revNewNodeIndex).end();
				assert(foundReverse);
				assert(revOldNodeOffset == 0);
				assert(revNewNodeOffset == params.graph.NodeLength(revNewNodeIndex)-1);
			}
		}
		assert(score == foundScore);
	}
#endif

};

#endif