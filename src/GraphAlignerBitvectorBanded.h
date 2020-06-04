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
#include "AlignmentCorrectnessEstimation.h"
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

	OnewayTrace getReverseTraceFromSeed(const std::string_view& sequence, int bigraphNodeId, size_t nodeOffset, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	{
		assert(!params.graph.HasPrefixSeeder());
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto initialBandwidth = BV::getInitialSliceExactPosition(params, bigraphNodeId, nodeOffset);
		auto slice = getSqrtSlices(sequence, initialBandwidth, numSlices, forceGlobal, reusableState);
		if (!params.preciseClipping && !forceGlobal) BV::removeWronglyAlignedEnd(slice);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}
		assert(sequence.size() <= std::numeric_limits<ScoreType>::max() - WordConfiguration<Word>::WordSize * 2);
		assert(slice.slices.back().minScore >= 0);
		assert(slice.slices.back().minScore <= (ScoreType)sequence.size() + (ScoreType)WordConfiguration<Word>::WordSize * 2);

		OnewayTrace result;
		if (params.preciseClipping)
		{
			result = BV::getReverseTraceFromTableExactEndPos(params, sequence, slice, reusableState);
		}
		else
		{
			result = BV::getReverseTraceFromTableStartLastRow(params, sequence, slice, reusableState);
		}

		return result;
	}

	OnewayTrace getBacktraceFullStart(const std::string_view& originalSequence, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	{
		assert(originalSequence.size() > 1);
		assert(!params.graph.HasPrefixSeeder());
		DPSlice startSlice;
		startSlice.j = -WordConfiguration<Word>::WordSize;
		startSlice.scores.addEmptyNodeMap(params.graph.NodeSize());
		startSlice.bandwidth = 1;
		startSlice.minScore = 0;
		startSlice.minScoreNode = 0;
		startSlice.minScoreNodeOffset = 0;
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
			node.endSlice = {0, 0, match ? 0 : 1};
			node.exists = true;
		}
		std::string_view alignableSequence { originalSequence.data()+1, originalSequence.size() - 1 };
		assert(alignableSequence.size() > 0);
		size_t numSlices = (alignableSequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		auto slice = getSqrtSlices(alignableSequence, startSlice, numSlices, forceGlobal, reusableState);
		if (!params.preciseClipping && !forceGlobal) BV::removeWronglyAlignedEnd(slice);
		if (slice.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}

		OnewayTrace result;
		if (params.preciseClipping)
		{
			result = BV::getReverseTraceFromTableExactEndPos(params, alignableSequence, slice, reusableState);
		}
		else
		{
			result = BV::getReverseTraceFromTableStartLastRow(params, alignableSequence, slice, reusableState);
		}
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

#ifdef EXTRACORRECTNESSASSERTIONS
	template <bool HasVectorMap, bool PreviousHasVectorMap>
	void checkNodeBoundaryCorrectness(const NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, const std::string_view& sequence, size_t j, ScoreType maxScore, ScoreType previousMaxScore) const
	{
		for (auto pair : currentSlice)
		{
			auto node = pair.first;
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
			for (auto neighbor : params.graph.inNeighbors[node])
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
				assert(pair.second.startSlice.getValue(0) == foundMinScore);
			}
			bool hasMatch = Common::characterMatch(sequence[0], params.graph.NodeSequences(node, 0));
			for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
			{
				if (j+i >= sequence.size()) break;
				eq = Common::characterMatch(sequence[j+i], params.graph.NodeSequences(node, 0));
				if (eq) hasMatch = true;
				ScoreType foundMinScore = pair.second.startSlice.getValue(i-1)+1;
				for (auto neighbor : params.graph.inNeighbors[node])
				{
					if (!currentSlice.hasNode(neighbor)) continue;
					if (!currentSlice.node(neighbor).exists) continue;
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i)+1);
					foundMinScore = std::min(foundMinScore, currentSlice.node(neighbor).endSlice.getValue(i-1) + (eq ? 0 : 1));
				}
				if (pair.second.startSlice.getValue(i) <= maxScore || foundMinScore <= maxScore)
				{
					assert(pair.second.startSlice.getValue(i) == foundMinScore);
				}
			}
		}
	}
#endif

	template <bool HasVectorMap, bool PreviousHasVectorMap, typename PriorityQueue>
	NodeCalculationResult calculateSlice(const std::string_view& sequence, size_t j, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& currentSlice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, std::vector<bool>& currentBand, const std::vector<bool>& previousBand, PriorityQueue& calculableQueue, ScoreType previousQuitScore, int bandwidth, ScoreType previousMinScore) const
	{
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

		assert(previousSlice.size() > 0);
		ScoreType zeroScore = previousMinScore*priorityMismatchPenalty - j - 64;
		if (j == 0)
		{
			for (auto node : previousSlice)
			{
				assert(node.second.minScore <= previousQuitScore);
				WordSlice startSlice = BV::getSourceSliceFromScore(node.second.startSlice.scoreEnd);
				if (calculableQueue.IsComponentPriorityQueue())
				{
					calculableQueue.insert(params.graph.componentNumber[node.first], node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
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
				assert(node.second.exists);
				if (node.second.minScore > previousQuitScore) continue;
				if (params.graph.linearizable[node.first])
				{
					auto neighbor = params.graph.inNeighbors[node.first][0];
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
					calculableQueue.insert(params.graph.componentNumber[node.first], node.second.minScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
				else
				{
					calculableQueue.insert(node.second.minScore*priorityMismatchPenalty - j - zeroScore, EdgeWithPriority { node.first, node.second.minScore - previousMinScore, startSlice, true });
				}
			}
		}
		assert(calculableQueue.size() > 0);
		
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
			if (!currentBand[i])
			{
				assert(!currentSlice.hasNode(i));
				currentSlice.addNode(i);
				currentBand[i] = true;
			}
			assert(currentBand[i]);
			const std::vector<EdgeWithPriority>* extras;
			extras = &calculableQueue.getExtras(i);
			auto& thisNode = currentSlice.node(i);
			auto oldEnd = thisNode.endSlice;
			if (!thisNode.exists) oldEnd = { 0, 0, std::numeric_limits<ScoreType>::max() };
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
			if (i < params.graph.firstAmbiguous)
			{
				if (params.preciseClipping)
				{
					nodeCalc = BV::calculateNodeClipPrecise(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.NodeChunks(i));
					assert(nodeCalc.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
				}
				else
				{
					nodeCalc = BV::calculateNodeClipApprox(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.NodeChunks(i));
				}
			}
			else
			{
				if (params.preciseClipping)
				{
					nodeCalc = BV::calculateNodeClipPrecise(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.AmbiguousNodeChunks(i));
					assert(nodeCalc.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
				}
				else
				{
					nodeCalc = BV::calculateNodeClipApprox(params, i, thisNode, EqV, previousThisNode, *extras, previousBand, params.graph.AmbiguousNodeChunks(i));
				}
			}
			calculableQueue.pop();
			if (!calculableQueue.IsComponentPriorityQueue())
			{
				calculableQueue.removeExtras(i);
			}
			assert(nodeCalc.minScore <= previousQuitScore + bandwidth + params.graph.SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize);
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

			if (newEnd.scoreEnd != oldEnd.scoreEnd || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				assert(newEndMinScore >= previousMinScore);
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				if (newEndMinScore <= currentMinScoreAtEndRow + bandwidth)
				{
					for (auto neighbor : params.graph.outNeighbors[i])
					{
						if (calculableQueue.IsComponentPriorityQueue())
						{
							calculableQueue.insert(params.graph.componentNumber[neighbor], newEndMinScore, EdgeWithPriority { neighbor, newEndMinScore - previousMinScore, newEnd, false });
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
			if (params.preciseClipping && nodeCalc.maxExactEndposScore > result.maxExactEndposScore)
			{
				result.maxExactEndposScore = nodeCalc.maxExactEndposScore;
				result.maxExactEndposNode = nodeCalc.maxExactEndposNode;
			}
			assert(result.minScore == currentMinScoreAtEndRow);
			result.cellsProcessed += nodeCalc.cellsProcessed;
			assert(nodeCalc.cellsProcessed > 0);
#ifdef SLICEVERBOSE
			result.nodesProcessed++;
#endif
			if (result.cellsProcessed > params.maxCellsPerSlice) break;
		}

#ifdef EXTRACORRECTNESSASSERTIONS
		checkNodeBoundaryCorrectness<HasVectorMap, PreviousHasVectorMap>(currentSlice, previousSlice, sequence, j, currentMinScoreAtEndRow + bandwidth, previousQuitScore);
#endif

		assert(result.minScoreNode != std::numeric_limits<LengthType>::max());

		if (!params.preciseClipping && j + WordConfiguration<Word>::WordSize > sequence.size())
		{
			BV::flattenLastSliceEnd(params, currentSlice, previousSlice, result, j, sequence);
		}

#ifdef SLICEVERBOSE
		std::cerr << "prefilternodes " << currentSlice.size() << " ";
#endif

		calculableQueue.clear();

		return result;
	}

	template <typename PriorityQueue>
	void fillDPSlice(const std::string_view& sequence, DPSlice& slice, const DPSlice& previousSlice, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, PriorityQueue& calculableQueue, int bandwidth) const
	{
		NodeCalculationResult sliceResult;
		assert((ScoreType)previousSlice.bandwidth < std::numeric_limits<ScoreType>::max());
		assert((ScoreType)previousSlice.bandwidth >= 0);
		assert(previousSlice.minScore < std::numeric_limits<ScoreType>::max() - (ScoreType)previousSlice.bandwidth);
		if (slice.scoresVectorMap.hasVectorMapCurrently())
		{
			if (previousSlice.scoresVectorMap.hasVectorMapCurrently())
			{
				sliceResult = calculateSlice<true, true>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scoresVectorMap, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
			}
			else
			{
				sliceResult = calculateSlice<true, false>(sequence, slice.j, slice.scoresVectorMap, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
			}
			slice.scores = slice.scoresVectorMap.getMapSlice();
		}
		else
		{
			assert(!previousSlice.scoresVectorMap.hasVectorMapCurrently());
			sliceResult = calculateSlice<false, false>(sequence, slice.j, slice.scores, previousSlice.scores, currentBand, previousBand, calculableQueue, previousSlice.minScore + previousSlice.bandwidth, bandwidth, previousSlice.minScore);
		}
		slice.cellsProcessed = sliceResult.cellsProcessed;
		slice.minScoreNode = sliceResult.minScoreNode;
		slice.minScoreNodeOffset = sliceResult.minScoreNodeOffset;
		slice.minScore = sliceResult.minScore;
		slice.maxExactEndposScore = sliceResult.maxExactEndposScore + slice.j;
		slice.maxExactEndposNode = sliceResult.maxExactEndposNode;
		assert(!params.preciseClipping || sliceResult.maxExactEndposScore != std::numeric_limits<ScoreType>::min());
		assert(!params.preciseClipping || sliceResult.maxExactEndposScore >= -((ScoreType)slice.j + WordConfiguration<Word>::WordSize) * 3);
		assert(!params.preciseClipping || sliceResult.maxExactEndposScore <= ((ScoreType)slice.j + WordConfiguration<Word>::WordSize) * 3);
		assert(!params.preciseClipping || slice.maxExactEndposScore <= ((ScoreType)slice.j + WordConfiguration<Word>::WordSize) * 3);
		assert(!params.preciseClipping || slice.maxExactEndposScore >= -((ScoreType)slice.j + WordConfiguration<Word>::WordSize) * 3);
		assert(slice.minScore >= previousSlice.minScore);
		slice.correctness = slice.correctness.NextState(slice.minScore - previousSlice.minScore, WordConfiguration<Word>::WordSize);
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
	DPSlice pickMethodAndExtendFill(const std::string_view& sequence, const DPSlice& previous, const std::vector<bool>& previousBand, std::vector<bool>& currentBand, std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem>& nodesliceMap, PriorityQueue& calculableQueue, int bandwidth) const
	{
		if (!params.lowMemory)
		{
			DPSlice bandTest { &nodesliceMap };
			bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
			bandTest.correctness = previous.correctness;
			fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth);
			return bandTest;
		}
		else
		{
			DPSlice bandTest;
			bandTest.scores.addEmptyNodeMap(previous.scores.size());
			bandTest.j = previous.j + WordConfiguration<Word>::WordSize;
			bandTest.correctness = previous.correctness;
			fillDPSlice(sequence, bandTest, previous, previousBand, currentBand, calculableQueue, bandwidth);
			return bandTest;
		}
	}

	DPTable getSqrtSlices(const std::string_view& sequence, const DPSlice& initialSlice, size_t numSlices, bool forceGlobal, AlignerGraphsizedState& reusableState) const
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
		assert(lastSlice.correctness.CurrentlyCorrect());
		DPSlice rampSlice = lastSlice;
		size_t rampRedoIndex = -1;
		size_t rampUntil = 0;
#ifndef NDEBUG
		volatile size_t debugLastProcessedSlice;
		// we want to keep this variable for debugging purposes
		// useless self-assignment to prevent unused variable compilation warning
		debugLastProcessedSlice = debugLastProcessedSlice;
#endif
		for (size_t slice = 0; slice < numSlices; slice++)
		{
			int bandwidth = (params.rampBandwidth > params.initialBandwidth && rampUntil >= slice) ? params.rampBandwidth : params.initialBandwidth;
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
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, (slice % 2 == 0) ? reusableState.evenNodesliceMap : reusableState.oddNodesliceMap, reusableState.componentQueue, bandwidth);
			}
			else
			{
				newSlice = pickMethodAndExtendFill(sequence, lastSlice, reusableState.previousBand, reusableState.currentBand, (slice % 2 == 0) ? reusableState.evenNodesliceMap : reusableState.oddNodesliceMap, reusableState.calculableQueue, bandwidth);
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

			if ((rampUntil == slice-1 || (rampUntil < slice && newSlice.correctness.CurrentlyCorrect() && newSlice.correctness.FalseFromCorrect())))
			{
				rampSlice = lastSlice.getMapSlice();
				rampRedoIndex = slice-1;
			}
			assert(newSlice.j == lastSlice.j + WordConfiguration<Word>::WordSize);

			cellsProcessed += newSlice.cellsProcessed;

			if (newSlice.cellsProcessed >= params.maxCellsPerSlice)
			{
				newSlice.scoresNotValid = true;
			}

			// if forcing global alignment, the whole read must be aligned -> don't allow breaking even if the alignment looks wrong
			if (!forceGlobal)
			{
				if (!newSlice.correctness.CorrectFromCorrect())
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
				if (!newSlice.correctness.CurrentlyCorrect() && rampUntil < slice && params.rampBandwidth > params.initialBandwidth)
				{
					for (auto node : newSlice.scores)
					{
						assert(reusableState.currentBand[node.first]);
						reusableState.currentBand[node.first] = false;
					}
					for (auto node : lastSlice.scores)
					{
						assert(reusableState.previousBand[node.first]);
						reusableState.previousBand[node.first] = false;
					}
					lastSlice.scoresVectorMap.removeVectorArray();
					newSlice.scoresVectorMap.removeVectorArray();
					rampUntil = slice;
					std::swap(slice, rampRedoIndex);
					std::swap(lastSlice, rampSlice);
					for (auto node : lastSlice.scores)
					{
						assert(!reusableState.previousBand[node.first]);
						reusableState.previousBand[node.first] = true;
					}
					if (slice == (size_t)-1)
					{
						result.slices.clear();
					}
					while (result.slices.size() > 1 && result.slices.back().j > slice * WordConfiguration<Word>::WordSize) result.slices.pop_back();
					assert(slice == (size_t)-1 || result.slices.size() == slice+2);
					assert(result.slices.back().j == lastSlice.j);
#ifdef SLICEVERBOSE
					std::cerr << " ramp to " << slice;
					std::cerr << std::endl;
					if (result.slices.size() > 0) std::cerr << " slices.back().j " << result.slices.back().j; else std::cerr << " slices.size() 0";
					std::cerr << std::endl;
#endif
					continue;
				}
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

};

#endif
