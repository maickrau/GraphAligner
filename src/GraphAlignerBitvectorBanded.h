#ifndef GraphAlignerBitvectorBanded_h
#define GraphAlignerBitvectorBanded_h

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <string_view>
#include "AlignmentGraph.h"
#include "NodeSlice.h"
#include "CommonUtils.h"
#include "GraphAlignerWrapper.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"
#include "GraphAlignerBitvectorCommon.h"
#include "GraphAlignerCommon.h"
#include "ArrayPriorityQueue.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorBanded
{
private:
	using BV = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	using WordSlice = typename BV::WordSlice;
	using EqVector = typename BV::EqVector;
	using EdgeWithPriority = typename Common::EdgeWithPriority;
	using DPSlice = typename BV::DPSlice;
	using DPTable = typename BV::DPTable;
	using NodeCalculationResult = typename BV::NodeCalculationResult;
	const Params& params;
public:

	GraphAlignerBitvectorBanded(const Params& params) :
	params(params)
	{
	}

	std::vector<OnewayTrace> getMultiseedTraces(const std::string_view& sequence, const std::string_view& bwSequence, const std::vector<ProcessedSeedHit>& seedHits, AlignerGraphsizedState& reusableState, std::vector<ScoreType>& sliceMaxScores) const
	{
		std::vector<OnewayTrace> traces = getMultiseedTracesOneWay(sequence, seedHits, reusableState, sliceMaxScores);
		std::vector<OnewayTrace> revTraces;
		for (const auto& fwTrace : traces)
		{
			auto reversePos = params.graph.GetReverseDigraphPosition(fwTrace.trace[0].DPposition.node, fwTrace.trace[0].DPposition.nodeOffset);
			assert(params.graph.GetReverseDigraphPosition(reversePos.first, reversePos.second).first == fwTrace.trace[0].DPposition.node);
			assert(params.graph.GetReverseDigraphPosition(reversePos.first, reversePos.second).second == fwTrace.trace[0].DPposition.nodeOffset);

			assert(fwTrace.trace[0].DPposition.seqPos < sequence.size());
			size_t seqPos = sequence.size() - 1 - fwTrace.trace[0].DPposition.seqPos;
			assert(seqPos < sequence.size());

			std::vector<ProcessedSeedHit> fakeBwSeeds;
			fakeBwSeeds.emplace_back(seqPos, reversePos.first);

			std::vector<ScoreType> fakeMaxScores;
			fakeMaxScores.resize(sliceMaxScores.size(), 0);
			auto bwTraces = getMultiseedTracesOneWay(bwSequence, fakeBwSeeds, reusableState, fakeMaxScores);
			for (size_t i = 0; i < bwTraces.size(); i++)
			{
				revTraces.emplace_back(std::move(bwTraces[i]));
			}
		}
		return revTraces;
	}

	std::vector<OnewayTrace> getMultiseedTracesOneWay(const std::string_view& sequence, const std::vector<ProcessedSeedHit>& seedHits, AlignerGraphsizedState& reusableState, std::vector<ScoreType>& sliceMaxScores) const
	{
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto initialSlice = BV::getInitialEmptySlice();
		auto slice = getMultiseedSlices(sequence, initialSlice, numSlices, reusableState, seedHits, sliceMaxScores);
		std::vector<OnewayTrace> results = BV::getLocalMaximaTracesFromTable(params, sequence, slice, reusableState, true, true, sliceMaxScores);
		removeDuplicateTraces(results);
		return results;
	}

	OnewayTrace getReverseTraceFromSeed(const std::string_view& sequence, int bigraphNodeId, size_t nodeOffset, int Xdropcutoff, AlignerGraphsizedState& reusableState) const
	{
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto alignmentBandwidth = BV::getInitialSliceExactPosition(params, bigraphNodeId, nodeOffset);
		auto slice = getSlices(sequence, alignmentBandwidth, numSlices, Xdropcutoff, reusableState);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}
		assert(sequence.size() <= std::numeric_limits<ScoreType>::max() - WordConfiguration<Word>::WordSize * 2);
		assert(slice.slices.back().minScore >= 0);
		assert(slice.slices.back().minScore <= (ScoreType)sequence.size() + (ScoreType)WordConfiguration<Word>::WordSize * 2);

		OnewayTrace result;
		result = BV::getReverseTraceFromTableExactEndPos(params, sequence, slice, reusableState, true, false);

		return result;
	}

	OnewayTrace getBacktraceFullStart(const std::string_view& originalSequence, int Xdropcutoff, AlignerGraphsizedState& reusableState) const
	{
		assert(originalSequence.size() > 1);
		DPSlice startSlice;
		startSlice.j = -WordConfiguration<Word>::WordSize;
		startSlice.scores.addEmptyNodeMap(params.graph.NodeSize());
		startSlice.bandwidth = 1;
		startSlice.minScore = 0;
		startSlice.minScoreNode = 0;
		startSlice.minScoreNodeOffset = 0;
		startSlice.maxExactEndposScore = -params.XscoreErrorCost;
		startSlice.maxExactEndposNode = 0;
		char firstChar = originalSequence[0];
		for (size_t i = 0; i < params.graph.NodeSize(); i++)
		{
			startSlice.scores.addNodeToMap(i);
			startSlice.scores.setMinScore(i, 0);
			auto& node = startSlice.scores.node(i);
			bool match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, 0));
			node.startSlice = {0, 0, match ? 0 : 1};
			node.minScore = match ? 0 : 1;
			for (size_t j = 1; j < params.graph.NodeLength(i); j++)
			{
				bool oldMatch = match;
				match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, j));
				if (oldMatch && !match)
				{
					node.HP[j / params.graph.SPLIT_NODE_SIZE] |= ((Word)1) << (j % params.graph.SPLIT_NODE_SIZE);
				}
				else if (match && !oldMatch)
				{
					node.HN[j / params.graph.SPLIT_NODE_SIZE] |= ((Word)1) << (j % params.graph.SPLIT_NODE_SIZE);
				}
				if (match) node.minScore = 0;
			}
			if (node.minScore == 0)
			{
				startSlice.maxExactEndposScore = 0;
				startSlice.maxExactEndposNode = i;
			}
			node.endSlice = {0, 0, match ? 0 : 1};
			node.exists = true;
		}
		std::string_view alignableSequence { originalSequence.data()+1, originalSequence.size() - 1 };
		assert(alignableSequence.size() > 0);
		size_t numSlices = (alignableSequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto slice = getSlices(alignableSequence, startSlice, numSlices, Xdropcutoff, reusableState);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}

		OnewayTrace result;
		result = BV::getReverseTraceFromTableExactEndPos(params, alignableSequence, slice, reusableState, true, false);
		for (size_t i = 0; i < result.trace.size(); i++)
		{
			result.trace[i].DPposition.seqPos += 1;
		}
		std::reverse(result.trace.begin(), result.trace.end());
		result.trace[0].sequenceCharacter = originalSequence[0];
		assert(result.trace[0].DPposition.seqPos == 0);
		return result;
	}

private:

	int traceAlignmentScore(const OnewayTrace& trace) const
	{
		assert(trace.trace[0].DPposition.seqPos >= trace.trace.back().DPposition.seqPos);
		ScoreType score = (trace.trace[0].DPposition.seqPos - trace.trace.back().DPposition.seqPos + 1)*100;
		score -= (ScoreType)trace.score * params.XscoreErrorCost;
		return score;
	}

	void removeDifferentTracePrefix(OnewayTrace& removeFrom, const OnewayTrace& removeWith) const
	{
		size_t lastEqualOffset = 0;
		assert(removeFrom.trace.size() > 0);
		assert(removeWith.trace.size() > 0);
		assert(removeFrom.trace.back() == removeWith.trace.back());
		while (lastEqualOffset+1 < removeFrom.trace.size() && lastEqualOffset+1 < removeWith.trace.size() && removeFrom.trace[removeFrom.trace.size()-1-(lastEqualOffset+1)] == removeWith.trace[removeWith.trace.size()-1-(lastEqualOffset+1)])
		{
			lastEqualOffset += 1;
		}
		assert(lastEqualOffset < removeFrom.trace.size());
		assert(lastEqualOffset > 0);
		removeFrom.trace.erase(removeFrom.trace.begin(), removeFrom.trace.end() - lastEqualOffset);
	}

	void removeDuplicateTraces(std::vector<OnewayTrace>& traces) const
	{
		if (traces.size() == 0) return;
		std::sort(traces.begin(), traces.end(), [this](const OnewayTrace& left, const OnewayTrace& right) {
			return traceAlignmentScore(left) > traceAlignmentScore(right);
		});
		for (size_t i = traces.size()-1; i > 0; i--)
		{
			for (size_t j = 0; j < i; j++)
			{
				if (traces[i].trace.back().DPposition == traces[j].trace.back().DPposition)
				{
					// times 100 because int scores are 1/100ths
					if (params.clipAmbiguousEnds >= 0 && traceAlignmentScore(traces[i]) >= traceAlignmentScore(traces[j]) - params.clipAmbiguousEnds * 100)
					{
						removeDifferentTracePrefix(traces[j], traces[i]);
					}
					std::swap(traces[i], traces.back());
					traces.pop_back();
					break;
				}
			}
		}
	}

#ifdef EXTRACORRECTNESSASSERTIONS
	template <bool HasVectorMap, bool PreviousHasVectorMap>
	void checkNodeBoundaryCorrectness(const NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, const std::string_view& sequence, size_t j, ScoreType maxScore, ScoreType previousMaxScore, const std::vector<bool>& hasSeedStart, const WordSlice seedstartSlice, const WordSlice fakeSlice) const
	{
		assert(previousMaxScore <= maxScore || seedstartSlice.getScoreBeforeStart() <= maxScore);
		for (auto pair : currentSlice)
		{
			auto node = pair.first;
			WordSlice extraSlice = hasSeedStart[node] ? seedstartSlice : fakeSlice;
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists)
			{
				assert(pair.second.startSlice.getScoreBeforeStart() <= previousSlice.node(node).startSlice.scoreEnd);
			}
			bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(node, 0));
			ScoreType foundMinScore = std::numeric_limits<ScoreType>::max();
			foundMinScore = std::min(foundMinScore, currentSlice.node(node).startSlice.getScoreBeforeStart()+1);
			if (previousSlice.hasNode(node) && previousSlice.node(node).exists && previousSlice.node(node).minScore <= previousMaxScore)
			{
				foundMinScore = std::min(foundMinScore, previousSlice.node(node).startSlice.scoreEnd+1);
			}
			for (auto neighbor : params.graph.InNeighbors(node))
			{
				if (currentSlice.hasNode(neighbor) && currentSlice.node(neighbor).exists)
				{
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(0)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getScoreBeforeStart() + (eq?0:1));
				}
				if (previousSlice.hasNode(neighbor) && previousSlice.node(neighbor).exists && previousSlice.node(neighbor).minScore <= previousMaxScore)
				{
					foundMinScore = std::min(foundMinScore, previousSlice.node(neighbor).endSlice.scoreEnd + (eq ? 0 : 1));
				}
			}
			if (pair.second.startSlice.getValue(0) <= maxScore || foundMinScore <= maxScore)
			{
				assert(foundMinScore != std::numeric_limits<ScoreType>::max());
				// todo figure out why breaks
				assert(pair.second.startSlice.getValue(0) == foundMinScore || pair.second.startSlice.getValue(0) == extraSlice.getValue(0));
			}
			bool hasMatch = Common::characterMatch(sequence[0], params.graph.NodeSequences(node, 0));
			for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
			{
				if (j+i >= sequence.size()) break;
				eq = Common::characterMatch(sequence[j+i], params.graph.NodeSequences(node, 0));
				if (eq) hasMatch = true;
				ScoreType foundMinScore = pair.second.startSlice.getValue(i-1)+1;
				for (auto neighbor : params.graph.InNeighbors(node))
				{
					if (!currentSlice.hasNode(neighbor)) continue;
					if (!currentSlice.node(neighbor).exists) continue;
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i-1) + (eq ? 0 : 1));
				}
				if (pair.second.startSlice.getValue(i) <= maxScore || foundMinScore <= maxScore)
				{
					// todo figure out why breaks
					assert(pair.second.startSlice.getValue(i) == foundMinScore || pair.second.startSlice.getValue(i) == extraSlice.getValue(i));
				}
			}
		}
	}
#endif

	template <bool HasVectorMap, bool PreviousHasVectorMap, typename PriorityQueue>
	void addSeedHitToScoresAndQueue(const ProcessedSeedHit& seedHit, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, PriorityQueue& calculableQueue, const WordSlice extraSlice, phmap::flat_hash_map<size_t, ScoreType>& nodeMaxExactEndposScore, bool storeNodeExactEndposScores) const
	{
#ifdef SLICEVERBOSE
		std::cerr << " " << seedHit.alignmentGraphNodeId << "(" << seedHit.nodeID << ")";
#endif
		size_t node = seedHit.alignmentGraphNodeId;
		assert(calculableQueue.IsComponentPriorityQueue());
		if (calculableQueue.IsComponentPriorityQueue())
		{
			EdgeWithPriority insertEdge { node, extraSlice.getValue(0), extraSlice, true };
			insertEdge.forceCalculation = true;
			calculableQueue.insert(params.graph.ComponentNumber(node), extraSlice.getValue(0), insertEdge);
		}
	}

	template <bool HasVectorMap, bool PreviousHasVectorMap, typename PriorityQueue>
	NodeCalculationResult calculateSlice(const std::string_view& sequence, const size_t j, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, PriorityQueue& calculableQueue, ScoreType previousQuitScore, ScoreType bandwidth, ScoreType previousMinScore, const std::vector<ProcessedSeedHit>& seedHits, size_t seedhitStart, size_t seedhitEnd, const WordSlice seedstartSlice, std::vector<bool>& hasSeedStart, std::unordered_set<size_t>& seedstartNodes, phmap::flat_hash_map<size_t, ScoreType>& nodeMaxExactEndposScore, bool storeNodeExactEndposScores, const std::vector<bool>& allowedBigraphNodes) const
	{
		if (previousMinScore == std::numeric_limits<ScoreType>::max() - bandwidth - 1)
		{
			assert(seedstartSlice.scoreEnd != std::numeric_limits<ScoreType>::max());
			previousMinScore = seedstartSlice.getScoreBeforeStart();
			previousQuitScore = seedstartSlice.getScoreBeforeStart() + bandwidth;
		}
		if (seedstartSlice.scoreEnd != std::numeric_limits<ScoreType>::max() && previousMinScore > seedstartSlice.scoreEnd)
		{
			previousMinScore = seedstartSlice.scoreEnd;
			previousQuitScore = std::max(previousQuitScore, seedstartSlice.scoreEnd + bandwidth);
		}
		double averageErrorRate = 0;
		if (j > 0)
		{
			averageErrorRate = (double)previousMinScore / (double)j;
		}
		double priorityMismatchPenalty = WordConfiguration<Word>::WordSize;
		if (averageErrorRate > 1.0/(double)WordConfiguration<Word>::WordSize)
		{
			priorityMismatchPenalty = 1.0 / averageErrorRate;
		}
		assert(priorityMismatchPenalty >= 0);
		assert(priorityMismatchPenalty <= WordConfiguration<Word>::WordSize);
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max() - bandwidth - 1;
		result.minScoreNode = std::numeric_limits<LengthType>::max();
		result.minScoreNodeOffset = std::numeric_limits<LengthType>::max();
		result.maxExactEndposNode = std::numeric_limits<LengthType>::min();
		result.maxExactEndposScore = std::numeric_limits<ScoreType>::min();
		result.cellsProcessed = 0;
#ifdef SLICEVERBOSE
		result.nodesProcessed = 0;
#endif

		EqVector EqV = BV::getEqVector(sequence, j);

		assert(previousSlice.size() > 0 || seedhitStart != std::numeric_limits<size_t>::max());
		ScoreType zeroScore = previousMinScore*priorityMismatchPenalty - j - 64;
		if (j == 0)
		{
			for (auto node : previousSlice)
			{
				assert(allowedBigraphNodes[params.graph.BigraphNodeID(node.first)]);
				assert(node.second.minScore <= previousQuitScore);
				WordSlice startSlice = BV::getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (calculableQueue.IsComponentPriorityQueue())
				{
					calculableQueue.insert(params.graph.ComponentNumber(node.first), node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore*priorityMismatchPenalty - j - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		else
		{
			for (auto node : previousSlice)
			{
				assert(allowedBigraphNodes[params.graph.BigraphNodeID(node.first)]);
				assert(node.second.exists);
				assert(node.second.minScore <= node.second.startSlice.scoreEnd);
				assert(node.second.minScore <= node.second.endSlice.scoreEnd);
				if (node.second.minScore > previousQuitScore) continue;
				if (params.graph.Linearizable(node.first))
				{
					auto neighbor = params.graph.InNeighbors(node.first)[0];
				 	if (previousBand[neighbor] && previousSlice.node(neighbor).endSlice.scoreEnd < previousQuitScore && previousSlice.node(neighbor).minScore < previousQuitScore)
				 	{
				 		//linear area, no need to add the later node into the queue 
				 		//because calculating the earlier node will guarantee that the later node will get added
				 		continue;
				 	}
				}
				WordSlice startSlice = BV::getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (calculableQueue.IsComponentPriorityQueue())
				{
					calculableQueue.insert(params.graph.ComponentNumber(node.first), node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore*priorityMismatchPenalty - j - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		assert(calculableQueue.size() > 0 || seedhitStart != std::numeric_limits<size_t>::max());

#ifdef SLICEVERBOSE
		std::cerr << "seeds";
#endif

		std::vector<size_t> clearSeedStarts;
		if (seedstartSlice.getScoreBeforeStart() < previousQuitScore + WordConfiguration<Word>::WordSize)
		{
#ifdef SLICEVERBOSE
			std::cerr << " " << seedhitEnd - seedhitStart << ":";
#endif
			for (size_t i = seedhitStart; i < seedhitEnd; i++)
			{
				if (hasSeedStart[seedHits[i].alignmentGraphNodeId]) continue;
				assert(allowedBigraphNodes[params.graph.BigraphNodeID(seedHits[i].alignmentGraphNodeId)]);
				addSeedHitToScoresAndQueue(seedHits[i], currentSlice, previousSlice, currentBand, previousBand, calculableQueue, seedstartSlice, nodeMaxExactEndposScore, storeNodeExactEndposScores);
				assert(!hasSeedStart[seedHits[i].alignmentGraphNodeId]);
				hasSeedStart[seedHits[i].alignmentGraphNodeId] = true;
				clearSeedStarts.push_back(seedHits[i].alignmentGraphNodeId);
				seedstartNodes.insert(seedHits[i].alignmentGraphNodeId);
			}
#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif
		}
		else
		{
#ifdef SLICEVERBOSE
			std::cerr << " not added because of score" << std::endl;
#endif
		}

		WordSlice fakeSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
		ScoreType currentMinScoreAtEndRow = result.minScore;
		while (calculableQueue.size() > 0)
		{
			auto pair = calculableQueue.top();
			if (!calculableQueue.IsComponentPriorityQueue())
			{
				if (pair.priority > currentMinScoreAtEndRow + bandwidth) break;
			}
			if (calculableQueue.extraSize(pair.target) == 0)
			{
				calculableQueue.pop();
				continue;
			}
			auto i = pair.target;
			bool firstCalc = false;
			assert(allowedBigraphNodes[params.graph.BigraphNodeID(i)]);
			if (!currentBand[i])
			{
				assert(!currentSlice.hasNode(i));
				currentSlice.addNode(i);
				currentBand[i] = true;
				firstCalc = true;
				if (hasSeedStart[i])
				{
					currentSlice.addNodeToMap(i);
					currentSlice.setMinScore(i, seedstartSlice.scoreEnd);
					auto& nodeScores = currentSlice.node(i);
					nodeScores.startSlice = seedstartSlice;
					nodeScores.endSlice = seedstartSlice;
					nodeScores.exists = true;
				}
			}
			assert(currentBand[i]);
			const std::vector<EdgeWithPriority>* extras;
			extras = &calculableQueue.getExtras(i);
			auto& thisNode = currentSlice.node(i);
			auto oldEnd = thisNode.endSlice;
			if (!thisNode.exists) oldEnd = { 0, 0, std::numeric_limits<ScoreType>::max() };
			WordSlice extraSlice = hasSeedStart[i] ? seedstartSlice : fakeSlice;
			if (extraSlice.scoreEnd != std::numeric_limits<ScoreType>::max()) oldEnd = oldEnd.mergeWith(extraSlice);
			typename NodeSlice<LengthType, ScoreType, Word, HasVectorMap>::NodeSliceMapItem previousThisNode;

			if (previousBand[i])
			{
				previousThisNode = previousSlice.node(i);
				assert(previousThisNode.exists);
			}
			else
			{
				for (size_t chunk = 0; chunk < previousThisNode.NUM_CHUNKS; chunk++)
				{
					previousThisNode.HP[chunk] = WordConfiguration<Word>::AllOnes;
					previousThisNode.HN[chunk] = WordConfiguration<Word>::AllZeros;
				}
				previousThisNode.exists = false;
			}
			NodeCalculationResult nodeCalc;
			if (i < params.graph.FirstAmbiguous())
			{
				nodeCalc = BV::calculateNodeClipPrecise(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.NodeChunks(i), extraSlice, j);
				assert(nodeCalc.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
			}
			else
			{
				nodeCalc = BV::calculateNodeClipPrecise(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.AmbiguousNodeChunks(i), extraSlice, j);
				assert(nodeCalc.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
			}
			assert(nodeCalc.minScore != std::numeric_limits<ScoreType>::max());
			calculableQueue.pop();
			if (!calculableQueue.IsComponentPriorityQueue())
			{
				calculableQueue.removeExtras(i);
			}
			if (storeNodeExactEndposScores)
			{
				nodeMaxExactEndposScore[i] = std::max(nodeCalc.maxExactEndposScore, nodeMaxExactEndposScore[i]);
			}
			// todo fix
			// assert(nodeCalc.minScore <= previousQuitScore + bandwidth + params.graph.SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize);
			assert(nodeCalc.minScore <= (ScoreType)j + bandwidth + (ScoreType)WordConfiguration<Word>::WordSize + (ScoreType)WordConfiguration<Word>::WordSize);
			currentMinScoreAtEndRow = std::min(currentMinScoreAtEndRow, nodeCalc.minScore);
			currentSlice.setMinScoreIfSmaller(i, nodeCalc.minScore);
#ifdef SLICEVERBOSE
			volatile size_t firstslices = currentSlice.node(i).firstSlicesCalcedWhenCalced;
			volatile size_t calcedslices = currentSlice.node(i).slicesCalcedWhenCalced;
			if (currentSlice.node(i).firstSlicesCalcedWhenCalced == std::numeric_limits<size_t>::max()) currentSlice.node(i).firstSlicesCalcedWhenCalced = result.cellsProcessed;
			if (currentSlice.node(i).slicesCalcedWhenCalced != std::numeric_limits<size_t>::max()) assert(currentSlice.node(i).slicesCalcedWhenCalced < result.cellsProcessed);
			currentSlice.node(i).slicesCalcedWhenCalced = result.cellsProcessed;
			assert(currentSlice.node(i).firstSlicesCalcedWhenCalced <= currentSlice.node(i).slicesCalcedWhenCalced);
#endif
			auto newEnd = thisNode.endSlice;

			if (firstCalc || newEnd.scoreEnd != oldEnd.scoreEnd || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				if (firstCalc) newEndMinScore = newEnd.getMinScore();
				// assert(newEndMinScore >= previousMinScore || newEndMinScore >= seedstartSlice.getScoreBeforeStart());
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				if (newEndMinScore <= currentMinScoreAtEndRow + bandwidth)
				{
					for (auto neighbor : params.graph.OutNeighbors(i))
					{
						if (!allowedBigraphNodes[params.graph.BigraphNodeID(neighbor)]) continue;
						if (calculableQueue.IsComponentPriorityQueue())
						{
							calculableQueue.insert(params.graph.ComponentNumber(neighbor), newEndMinScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
						}
						else
						{
							ScoreType newEndPriorityScore = newEnd.getChangedPriorityScore(oldEnd, j, priorityMismatchPenalty);
							assert(newEndPriorityScore != std::numeric_limits<ScoreType>::max());
							assert(newEndPriorityScore >= zeroScore);
							calculableQueue.insert(newEndPriorityScore - zeroScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
						}
					}
				}
			}
			if (nodeCalc.minScore < result.minScore)
			{
				result.minScore = nodeCalc.minScore;
				result.minScoreNode = nodeCalc.minScoreNode;
				result.minScoreNodeOffset = nodeCalc.minScoreNodeOffset;
			}
			if (nodeCalc.maxExactEndposScore > result.maxExactEndposScore)
			{
				result.maxExactEndposScore = nodeCalc.maxExactEndposScore;
				result.maxExactEndposNode = nodeCalc.maxExactEndposNode;
			}
			assert(currentSlice.node(i).minScore <= currentSlice.node(i).startSlice.scoreEnd);
			assert(currentSlice.node(i).minScore <= currentSlice.node(i).endSlice.scoreEnd);
			assert(result.minScore == currentMinScoreAtEndRow);
			result.cellsProcessed += nodeCalc.cellsProcessed;
			assert(nodeCalc.cellsProcessed > 0);
#ifdef SLICEVERBOSE
			result.nodesProcessed++;
#endif
			if (result.cellsProcessed > params.maxCellsPerSlice) break;
		}

#ifdef EXTRACORRECTNESSASSERTIONS
		checkNodeBoundaryCorrectness<HasVectorMap, PreviousHasVectorMap>(currentSlice, previousSlice, sequence, j, currentMinScoreAtEndRow + bandwidth, previousQuitScore, hasSeedStart, seedstartSlice, fakeSlice);
#endif
		for (auto node : currentSlice)
		{
			assert(node.second.exists);
			assert(node.second.minScore <= node.second.startSlice.scoreEnd);
			assert(node.second.minScore <= node.second.endSlice.scoreEnd);
		}

		for (auto node : clearSeedStarts)
		{
			assert(hasSeedStart[node]);
			hasSeedStart[node] = false;
		}

		assert(result.minScoreNode != std::numeric_limits<LengthType>::max() || (seedhitStart == seedhitEnd && seedhitStart != std::numeric_limits<size_t>::max()));

#ifdef SLICEVERBOSE
		std::cerr << "prefilternodes " << currentSlice.size() << " ";
#endif

		calculableQueue.clear();

		return result;
	}

	template <typename PriorityQueue>
	void fillDPSlice(const std::string_view& sequence, DPSlice& slice, const DPSlice& previousSlice, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, PriorityQueue& calculableQueue, ScoreType bandwidth, const std::vector<ProcessedSeedHit>& seedHits, size_t seedhitStart, size_t seedhitEnd, const WordSlice extraSlice, std::vector<bool>& hasSeedStart, bool storeNodeExactEndposScores, const std::vector<bool>& allowedBigraphNodes) const
	{
		NodeCalculationResult sliceResult;
		assert((ScoreType)previousSlice.bandwidth < std::numeric_limits<ScoreType>::max());
		assert((ScoreType)previousSlice.bandwidth >= 0);
		assert(previousSlice.minScore < std::numeric_limits<ScoreType>::max() - (ScoreType)previousSlice.bandwidth);
		if (slice.scoresVectorMap.hasVectorMapCurrently())
		{
			if (previousSlice.scoresVectorMap.hasVectorMapCurrently())
			{
				sliceResult = calculateSlice<true, true>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scoresVectorMap, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, slice.seedstartNodes, slice.nodeMaxExactEndposScore, storeNodeExactEndposScores, allowedBigraphNodes);
			}
			else
			{
				sliceResult = calculateSlice<true, false>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, slice.seedstartNodes, slice.nodeMaxExactEndposScore, storeNodeExactEndposScores, allowedBigraphNodes);
			}
			slice.scores = slice.scoresVectorMap.getMapSlice();
		}
		else
		{
			assert(!previousSlice.scoresVectorMap.hasVectorMapCurrently());
			sliceResult = calculateSlice<false, false>(sequence, slice.j, slice.scores, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, slice.seedstartNodes, slice.nodeMaxExactEndposScore, storeNodeExactEndposScores, allowedBigraphNodes);
		}
		slice.cellsProcessed = sliceResult.cellsProcessed;
		slice.minScoreNode = sliceResult.minScoreNode;
		slice.minScoreNodeOffset = sliceResult.minScoreNodeOffset;
		slice.minScore = sliceResult.minScore;
		slice.maxExactEndposScore = sliceResult.maxExactEndposScore;
		slice.maxExactEndposNode = sliceResult.maxExactEndposNode;
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || sliceResult.minScore <= (ScoreType)slice.j + (ScoreType)WordConfiguration<Word>::WordSize);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || sliceResult.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)sliceResult.maxExactEndposScore >= -((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)sliceResult.maxExactEndposScore <= ((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)slice.maxExactEndposScore <= ((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert((seedHits.size() != 0 && seedhitEnd == seedhitStart) || (double)slice.maxExactEndposScore >= -((double)slice.j + WordConfiguration<Word>::WordSize) * params.XscoreErrorCost);
		assert(slice.minScore >= previousSlice.minScore || slice.minScore >= extraSlice.getValue(0) || previousSlice.minScore == std::numeric_limits<ScoreType>::max() - bandwidth - 1);
		slice.bandwidth = bandwidth;
#ifdef SLICEVERBOSE
		slice.nodesProcessed = sliceResult.nodesProcessed;
		for (auto node : slice.scores)
		{
			if (currentBand[node.first])
			{
				slice.numCells += params.graph.NodeLength(node.first);
			}
		}
#endif
	}

	template <typename PriorityQueue>
	DPSlice pickMethodAndExtendFill(const std::string_view& sequence, const DPSlice& previous, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, PriorityQueue& calculableQueue, ScoreType bandwidth, const std::vector<ProcessedSeedHit>& seedHits, size_t seedhitStart, size_t seedhitEnd, const WordSlice extraSlice, std::vector<bool>& hasSeedStart, bool storeNodeExactEndposScores, const std::vector<bool>& allowedBigraphNodes) const
	{
		DPSlice bandTest;
		bandTest.scores.addEmptyNodeMap(previous.scores.size());
		bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
		fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth, seedHits, seedhitStart, seedhitEnd, extraSlice, hasSeedStart, storeNodeExactEndposScores, allowedBigraphNodes);
		return bandTest;
	}

	DPTable getSlices(const std::string_view& sequence, const DPSlice& initialSlice, size_t numSlices, int Xdropcutoff, AlignerGraphsizedState& reusableState) const
	{
		return getXdropSlices(sequence, initialSlice, numSlices, Xdropcutoff, reusableState);
	}

	DPTable getXdropSlices(const std::string_view& sequence, const DPSlice& initialSlice, size_t numSlices, double Xdropcutoff, AlignerGraphsizedState& reusableState) const
	{
		assert(initialSlice.j == (size_t)-WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		DPTable result;
		result.slices.reserve(numSlices + 1);
		size_t cellsProcessed = 0;
		std::vector<size_t> partOfComponent;
		{
			for (auto node : initialSlice.scores)
			{
				reusableState.previousBand[node.first] = true;
			}
		}
#ifndef NDEBUG
		debugLastRowMinScore = 0;
#endif
		DPSlice lastSlice = initialSlice;
		result.slices.push_back(initialSlice);
		ScoreType bestXScore = initialSlice.maxExactEndposScore;
		assert(bestXScore != std::numeric_limits<ScoreType>::min());
#ifndef NDEBUG
		volatile size_t debugLastProcessedSlice;
		// we want to keep this variable for debugging purposes
		// useless self-assignment to prevent unused variable compilation warning
		debugLastProcessedSlice = debugLastProcessedSlice;
#endif
		std::vector<ProcessedSeedHit> fakeSeeds;
		WordSlice fakeSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = params.alignmentBandwidth;
#ifndef NDEBUG
			debugLastProcessedSlice = slice;
			debugLastRowMinScore = lastSlice.minScore;
#endif
#ifdef SLICEVERBOSE
			auto timeStart = std::chrono::system_clock::now();
#endif
			DPSlice newSlice;
			if (reusableState.componentQueue.valid())
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.componentQueue, bandwidth, fakeSeeds, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), fakeSlice, reusableState.hasSeedStart, false, reusableState.allowedBigraphNodes);
			}
			else
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.calculableQueue, bandwidth, fakeSeeds, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), fakeSlice, reusableState.hasSeedStart, false, reusableState.allowedBigraphNodes);
			}
#ifdef SLICEVERBOSE
			auto timeEnd = std::chrono::system_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
			std::cerr << "slice " << slice << " bandwidth " << bandwidth << " minscore " << newSlice.minScore << " diff " << (newSlice.minScore - lastSlice.minScore) << " time " << time << " nodes " << newSlice.scores.size() << " slices " << newSlice.numCells << " nodesprocessed " << newSlice.nodesProcessed << " cellsprocessed " << newSlice.cellsProcessed << " overhead " << (100 * (int)(newSlice.cellsProcessed - newSlice.numCells) / (int)(newSlice.numCells)) << "%";
#endif
			assert(newSlice.minScore != std::numeric_limits<ScoreType>::max());
			assert(newSlice.minScoreNode != std::numeric_limits<LengthType>::max());
			assert(newSlice.minScoreNodeOffset != std::numeric_limits<LengthType>::max());
			assert(newSlice.scores.hasNode(newSlice.minScoreNode));
			assert(newSlice.minScoreNodeOffset < params.graph.NodeLength(newSlice.minScoreNode));
			assert(newSlice.maxExactEndposScore != std::numeric_limits<ScoreType>::min());

			if (newSlice.maxExactEndposScore > bestXScore) bestXScore = newSlice.maxExactEndposScore;

			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed >= params.maxCellsPerSlice)
			{
				newSlice.scoresNotValid = true;
			}

			if (newSlice.maxExactEndposScore < bestXScore - Xdropcutoff)
			{
#ifndef NDEBUG
				debugLastProcessedSlice = slice-1;
#endif
				for (auto node : lastSlice.scores)
				{
					assert(reusableState.previousBand[node.first]);
					reusableState.previousBand[node.first] = false;
				}
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
				lastSlice.scoresVectorMap.removeVectorArray();
				newSlice.scoresVectorMap.removeVectorArray();
				break;
			}

#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			result.slices.push_back(newSlice.getMapSlice());
			for (auto node : lastSlice.scores)
			{
				assert(reusableState.previousBand[node.first]);
				reusableState.previousBand[node.first] = false;
			}
			assert(newSlice.minScore != std::numeric_limits<ScoreType>::max());
			assert(newSlice.minScore >= lastSlice.minScore);
			if (slice == numSlices - 1)
			{
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice.scoresVectorMap.removeVectorArray();
			lastSlice = std::move(newSlice);
		}
		lastSlice.scoresVectorMap.removeVectorArray();

		assert(result.slices.size() <= numSlices + 1);

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
		}
#endif

#ifndef NDEBUG
		if (result.slices.size() > 0)
		{
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].j == result.slices[i-1].j + WordConfiguration<Word>::WordSize);
			}
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].minScore >= result.slices[i-1].minScore);
			}
		}
#endif
		return result;
	}

	DPTable getMultiseedSlices(const std::string_view& sequence, const DPSlice& initialSlice, size_t numSlices, AlignerGraphsizedState& reusableState, const std::vector<ProcessedSeedHit>& seedHits, const std::vector<ScoreType>& sliceMaxScores) const
	{
		assert(reusableState.componentQueue.valid());
		assert(initialSlice.j == (size_t)-WordConfiguration<Word>::WordSize);
		assert(initialSlice.j + numSlices * WordConfiguration<Word>::WordSize <= sequence.size() + WordConfiguration<Word>::WordSize);
		assert(numSlices >= 1);
		assert(numSlices == 1 || initialSlice.j + (numSlices-1) * WordConfiguration<Word>::WordSize <= sequence.size());
		DPTable result;
		result.slices.reserve(numSlices + 1);
		size_t cellsProcessed = 0;
		DPSlice lastSlice = initialSlice;
		result.slices.push_back(initialSlice);
		size_t lastSeedHit = 0;
		ScoreType XDropCurrentBest = 0;
		assert(params.Xdropcutoff > 0);
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = params.alignmentBandwidth;
#ifdef SLICEVERBOSE
			auto timeStart = std::chrono::system_clock::now();
#endif
			size_t nextSeedHit = lastSeedHit;
			assert(nextSeedHit == seedHits.size() || seedHits[nextSeedHit].seqPos / WordConfiguration<Word>::WordSize >= (lastSlice.j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize);
			while (nextSeedHit < seedHits.size() && seedHits[nextSeedHit].seqPos / WordConfiguration<Word>::WordSize == (lastSlice.j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize)
			{
				nextSeedHit += 1;
			}
			assert(nextSeedHit == seedHits.size() || seedHits[nextSeedHit].seqPos / WordConfiguration<Word>::WordSize > (lastSlice.j + WordConfiguration<Word>::WordSize) / WordConfiguration<Word>::WordSize);
			size_t seqOffset = lastSlice.j + WordConfiguration<Word>::WordSize;
			WordSlice seedSlice = BV::getSeedSlice(seqOffset, sequence.size(), params);
			assert(seedSlice.maxXScore(seqOffset, params.XscoreErrorCost) >= -(ScoreType)WordConfiguration<Word>::WordSize*100);
			assert(seedSlice.maxXScore(seqOffset, params.XscoreErrorCost) <= (ScoreType)WordConfiguration<Word>::WordSize*100);
			ScoreType possibleScoreRemaining = sequence.size() - seqOffset;
			assert(possibleScoreRemaining > 0);
			// // can't get an alignment which would be backtraced. fake set the seed set to be empty
			// if (possibleScoreRemaining < sliceMaxScores[slice]) lastSeedHit = nextSeedHit;
			DPSlice newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, reusableState.componentQueue, bandwidth, seedHits, lastSeedHit, nextSeedHit, seedSlice, reusableState.hasSeedStart, true, reusableState.allowedBigraphNodes);
			lastSeedHit = nextSeedHit;
			XDropCurrentBest = std::max(XDropCurrentBest, newSlice.maxExactEndposScore);
			if (newSlice.maxExactEndposScore < XDropCurrentBest - params.Xdropcutoff)
			{
				for (auto node : newSlice.scores)
				{
					reusableState.currentBand[node.first] = false;
				}
				newSlice.minScore = seedSlice.getScoreBeforeStart();
				newSlice.minScoreNode = std::numeric_limits<LengthType>::max();
				newSlice.maxExactEndposScore = std::numeric_limits<ScoreType>::min();
				newSlice.maxExactEndposNode = std::numeric_limits<LengthType>::max();
				newSlice.nodeMaxExactEndposScore.clear();
				newSlice.seedstartNodes.clear();
				newSlice.scores.clear();
				XDropCurrentBest = 0;
			}
#ifdef SLICEVERBOSE
			auto timeEnd = std::chrono::system_clock::now();
			auto time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
			std::cerr << "slice " << slice << " bandwidth " << bandwidth << " minscore " << newSlice.minScore << " diff " << (newSlice.minScore - lastSlice.minScore) << " maxX " << newSlice.maxExactEndposScore << " time " << time << " nodes " << newSlice.scores.size() << " slices " << newSlice.numCells << " nodesprocessed " << newSlice.nodesProcessed << " cellsprocessed " << newSlice.cellsProcessed << " overhead " << (100 * (int)(newSlice.cellsProcessed - newSlice.numCells) / (std::max((int)newSlice.numCells, (int)1))) << "%";
#endif

			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed >= params.maxCellsPerSlice)
			{
				newSlice.scoresNotValid = true;
			}

#ifdef SLICEVERBOSE
			std::cerr << std::endl;
#endif

			result.slices.push_back(newSlice.getMapSlice());
			for (auto node : lastSlice.scores)
			{
				assert(reusableState.previousBand[node.first]);
				reusableState.previousBand[node.first] = false;
			}
			if (slice == numSlices - 1)
			{
				for (auto node : newSlice.scores)
				{
					assert(reusableState.currentBand[node.first]);
					reusableState.currentBand[node.first] = false;
				}
			}
			else
			{
				std::swap(reusableState.previousBand, reusableState.currentBand);
			}
			lastSlice.scoresVectorMap.removeVectorArray();
			lastSlice = std::move(newSlice);
		}
		assert(lastSeedHit == seedHits.size());
		lastSlice.scoresVectorMap.removeVectorArray();

		assert(result.slices.size() == numSlices + 1);

#ifdef EXTRACORRECTNESSASSERTIONS
		assert(reusableState.calculableQueue.size() == 0);
		for (size_t i = 0; i < reusableState.currentBand.size(); i++)
		{
			assert(!reusableState.currentBand[i]);
			assert(!reusableState.previousBand[i]);
		}
#endif

#ifndef NDEBUG
		if (result.slices.size() > 0)
		{
			for (size_t i = 1; i < result.slices.size(); i++)
			{
				assert(result.slices[i].j == result.slices[i-1].j + WordConfiguration<Word>::WordSize);
			}
		}
#endif
		return result;
	}

};

#endif
