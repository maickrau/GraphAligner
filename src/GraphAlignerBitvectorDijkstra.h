#ifndef GraphAlignerBitvectorDijkstra_h
#define GraphAlignerBitvectorDijkstra_h

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorDijkstra
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

	GraphAlignerBitvectorDijkstra(const Params& params) :
	params(params)
	{
	}

	// OnewayTrace getReverseTraceFromSeed(const std::string_view& sequence, int bigraphNodeId, size_t nodeOffset, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	// {
	// 	size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
	// 	auto initialBandwidth = BV::getInitialSliceExactPosition(params, bigraphNodeId, nodeOffset);
	// 	auto slice = getTable(sequence, initialBandwidth, reusableState);
	// 	if (!params.preciseClipping && !forceGlobal) BV::removeWronglyAlignedEnd(slice);
	// 	if (slice.slices.size() <= 1)
	// 	{
	// 		return OnewayTrace::TraceFailed();
	// 	}
	// 	assert(sequence.size() <= std::numeric_limits<ScoreType>::max() - WordConfiguration<Word>::WordSize * 2);
	// 	assert(slice.slices.back().minScore >= 0);
	// 	assert(slice.slices.back().minScore <= (ScoreType)sequence.size() + (ScoreType)WordConfiguration<Word>::WordSize * 2);

	// 	OnewayTrace result;
	// 	if (params.preciseClipping)
	// 	{
	// 		result = BV::getReverseTraceFromTableExactEndPos(params, sequence, slice, reusableState);
	// 	}
	// 	else
	// 	{
	// 		result = BV::getReverseTraceFromTableStartLastRow(params, sequence, slice, reusableState);
	// 	}

	// 	return result;
	// }

	OnewayTrace getBacktraceFullStart(const std::string_view& originalSequence, bool forceGlobal, AlignerGraphsizedState& reusableState) const
	{
		assert(originalSequence.size() > 1);
		std::string_view alignableSequence { originalSequence.data()+1, originalSequence.size() - 1 };
		assert(alignableSequence.size() > 0);
		size_t numSlices = (alignableSequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		DPTable table;
		table.slices.resize(numSlices + 1);
		table.slices[0].j = -WordConfiguration<Word>::WordSize;
		table.slices[0].scores.addEmptyNodeMap(params.graph.NodeSize());
		table.slices[0].bandwidth = 1;
		table.slices[0].minScore = 1;
		table.slices[0].minScoreNode = 0;
		table.slices[0].minScoreNodeOffset = 0;
		char firstChar = originalSequence[0];
		for (size_t i = 0; i < params.graph.NodeSize(); i++)
		{
			table.slices[0].scores.addNodeToMap(i);
			table.slices[0].scores.setMinScore(i, 0);
			auto& node = table.slices[0].scores.node(i);
			bool match = Common::characterMatch(firstChar, params.graph.NodeSequences(i, 0));
			node.startSlice = {0, 0, match ? 0 : 1};
			assert(table.slices[0].scores.node(i).startSlice.scoreEnd == (match ? 0 : 1));
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
				table.slices[0].minScore = 0;
				table.slices[0].minScoreNode = i;
				// todo fix, but probably doesn't matter in practice
				table.slices[0].minScoreNodeOffset = 0;
			}
			node.endSlice = {0, 0, match ? 0 : 1};
			assert(table.slices[0].scores.node(i).endSlice.scoreEnd == (match ? 0 : 1));
			node.exists = true;
		}
		for (size_t i = 1; i < table.slices.size(); i++)
		{
			table.slices[i].j = (i - 1) * WordConfiguration<Word>::WordSize;
			table.slices[i].scores.addEmptyNodeMap(params.graph.NodeSize());
			table.slices[i].bandwidth = 1;
			table.slices[i].minScore = table.slices[i].j + WordConfiguration<Word>::WordSize;
			table.slices[i].minScoreNode = 0;
			table.slices[i].minScoreNodeOffset = 0;
		}
		fillTable(table, alignableSequence, reusableState);
		if (!params.preciseClipping && !forceGlobal) BV::removeWronglyAlignedEnd(table);
		if (table.slices.size() <= 1)
		{
			return OnewayTrace::TraceFailed();
		}

		OnewayTrace result;
		if (params.preciseClipping)
		{
			result = BV::getReverseTraceFromTableExactEndPos(params, alignableSequence, table, reusableState);
		}
		else
		{
			result = BV::getReverseTraceFromTableStartLastRow(params, alignableSequence, table, reusableState);
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

	void fillTable(DPTable& table, const std::string_view& sequence, AlignerGraphsizedState& reusableState) const
	{
		assert(!params.preciseClipping);
		assert(reusableState.dijkstraQueue.size() == 0);
		assert(table.slices.size() > 0);
		for (auto node : table.slices[0].scores)
		{
			WordSlice startSlice = BV::getSourceSliceFromScore(node.second.startSlice.scoreEnd);
			WordSlice endSlice = BV::getSourceSliceFromScore(node.second.endSlice.scoreEnd);
			assert(node.second.startSlice.scoreEnd == 0 || node.second.startSlice.scoreEnd == 1);
			reusableState.dijkstraQueue.insert(table.slices[0].scores.minScore(node.first), EdgeWithPriority { node.first, table.slices[0].scores.minScore(node.first), startSlice, true, 0 });
			for (auto neighbor : params.graph.outNeighbors[node.first])
			{
				assert(node.second.endSlice.scoreEnd == 0 || node.second.endSlice.scoreEnd == 1);
				reusableState.dijkstraQueue.insert(node.second.endSlice.scoreEnd, EdgeWithPriority { neighbor, node.second.endSlice.scoreEnd, endSlice, false, 0 });
			}
		}
		size_t numSlices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		assert(table.slices.size() == numSlices+1);
		std::vector<EqVector> EqV;
		size_t zeroScore = 0;
		for (size_t i = 0; i < numSlices; i++)
		{
			EqV.emplace_back(BV::getEqVector(sequence, i * WordConfiguration<Word>::WordSize));
		}
		size_t lastRowScore = std::numeric_limits<size_t>::max();
		assert(reusableState.dijkstraQueue.zero() == 0);
		while (true)
		{
			if (reusableState.dijkstraQueue.size() == 0) break;
			auto edge = reusableState.dijkstraQueue.top();
			if (reusableState.dijkstraQueue.extraSize(edge.slice, edge.target) == 0)
			{
				reusableState.dijkstraQueue.pop();
				continue;
			}
			if (edge.priority != zeroScore)
			{
				assert(edge.priority > zeroScore);
				reusableState.dijkstraQueue.increaseScore(edge.priority - zeroScore);
				zeroScore = edge.priority;
			}
			if (zeroScore >= lastRowScore) break;
			size_t tableSlice = edge.slice+1;
			assert(tableSlice > 0);
			assert(tableSlice <= numSlices);
			size_t slice = tableSlice-1;
			const std::vector<EdgeWithPriority>* extras;
			size_t i = edge.target;
			extras = &reusableState.dijkstraQueue.getExtras(slice, i);
			assert(extras->size() > 0);
			if (!table.slices[tableSlice].scores.hasNode(i))
			{
				table.slices[tableSlice].scores.addNodeToMap(i);
				table.slices[tableSlice].scores.setMinScore(i, zeroScore + 64);
			}
			auto& thisNode = table.slices[tableSlice].scores.node(i);
			auto oldEnd = thisNode.endSlice;
			auto oldHP = thisNode.HP[0];
			auto oldHN = thisNode.HN[0];
			if (!thisNode.exists) oldEnd = { 0, 0, std::numeric_limits<ScoreType>::max() };
			typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem previousThisNode;

			if (table.slices[tableSlice-1].scores.hasNode(i))
			{
				previousThisNode = table.slices[tableSlice-1].scores.node(i);
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
				nodeCalc = BV::calculateNodeClipApprox(params, i, thisNode, EqV[slice], previousThisNode, *extras, table.slices[tableSlice-1], params.graph.NodeChunks(i));
			}
			else
			{
				nodeCalc = BV::calculateNodeClipApprox(params, i, thisNode, EqV[slice], previousThisNode, *extras, table.slices[tableSlice-1], params.graph.AmbiguousNodeChunks(i));
			}
			if (tableSlice == numSlices)
			{
				if (nodeCalc.minScore < lastRowScore)
				{
					lastRowScore = nodeCalc.minScore;
				}
			}
			if (nodeCalc.minScore < table.slices[tableSlice].minScore)
			{
				table.slices[tableSlice].minScore = nodeCalc.minScore;
				table.slices[tableSlice].minScoreNode = nodeCalc.minScoreNode;
				table.slices[tableSlice].minScoreNodeOffset = nodeCalc.minScoreNodeOffset;
			}
			reusableState.dijkstraQueue.pop();
			reusableState.dijkstraQueue.removeExtras(slice, i);
			assert(nodeCalc.minScore <= zeroScore + params.graph.SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize);
			table.slices[tableSlice].scores.setMinScoreIfSmaller(i, nodeCalc.minScore);
#ifdef SLICEVERBOSE
			volatile size_t firstslices = table.slices[tableSlice].node(i).firstSlicesCalcedWhenCalced;
			volatile size_t calcedslices = table.slices[tableSlice].node(i).slicesCalcedWhenCalced;
			if (table.slices[tableSlice].node(i).firstSlicesCalcedWhenCalced == std::numeric_limits<size_t>::max()) table.slices[tableSlice].node(i).firstSlicesCalcedWhenCalced = result.cellsProcessed;
			if (table.slices[tableSlice].node(i).slicesCalcedWhenCalced != std::numeric_limits<size_t>::max()) assert(table.slices[tableSlice].node(i).slicesCalcedWhenCalced < result.cellsProcessed);
			table.slices[tableSlice].node(i).slicesCalcedWhenCalced = table.cellsProcessed;
			assert(table.slices[tableSlice].node(i).firstSlicesCalcedWhenCalced <= table.slices[tableSlice].node(i).slicesCalcedWhenCalced);
#endif
			auto newEnd = thisNode.endSlice;
			auto newHP = thisNode.HP[0];
			auto newHN = thisNode.HN[0];

			if (newEnd.scoreEnd != oldEnd.scoreEnd || newHP != oldHP || newHN != oldHN)
			{
				ScoreType newEndMinScore = changedHorizontal(newEnd, newHP, newHN, oldEnd, oldHP, oldHN, params.graph.NodeLength(i));
				assert(newEndMinScore >= zeroScore);
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				reusableState.dijkstraQueue.insert(newEndMinScore, EdgeWithPriority { i, newEndMinScore, BV::getSourceSliceFromScore(thisNode.startSlice.scoreEnd), true, slice+1 });
			}
			if (newEnd.scoreEnd != oldEnd.scoreEnd)
			{
				ScoreType newEndMinScore = newEnd.scoreEnd;
				assert(newEndMinScore >= zeroScore);
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				for (auto neighbor : params.graph.outNeighbors[i])
				{
					reusableState.dijkstraQueue.insert(newEndMinScore, EdgeWithPriority { neighbor, newEndMinScore, BV::getSourceSliceFromScore(newEnd.scoreEnd), false, slice+1 });
				}
			}
			if (newEnd.scoreEnd != oldEnd.scoreEnd || newEnd.VP != oldEnd.VP || newEnd.VN != oldEnd.VN)
			{
				ScoreType newEndMinScore = newEnd.changedMinScore(oldEnd);
				assert(newEndMinScore >= zeroScore);
				assert(newEndMinScore != std::numeric_limits<ScoreType>::max());
				for (auto neighbor : params.graph.outNeighbors[i])
				{
					reusableState.dijkstraQueue.insert(newEndMinScore, EdgeWithPriority { neighbor, newEndMinScore, newEnd, false, slice });
				}
			}
			table.slices[tableSlice].cellsProcessed += nodeCalc.cellsProcessed;
			assert(nodeCalc.cellsProcessed > 0);
#ifdef SLICEVERBOSE
			table.slices[tableSlice].nodesProcessed++;
#endif
		}
		reusableState.dijkstraQueue.clear();
	}

	ScoreType changedHorizontal(WordSlice newEnd, Word newHP, Word newHN, WordSlice oldEnd, Word oldHP, Word oldHN, size_t size) const
	{
		ScoreType newScore = newEnd.scoreEnd;
		ScoreType oldScore = oldEnd.scoreEnd;
		ScoreType result = std::numeric_limits<ScoreType>::max();
		if (newScore < oldScore) result = newScore;
		for (size_t i = size - 1; i > 0; i--)
		{
			newScore += (newHN >> i) & 1;
			newScore -= (newHP >> i) & 1;
			oldScore += (oldHN >> i) & 1;
			oldScore -= (oldHP >> i) & 1;
			if (newScore < oldScore) result = std::min(result, newScore);
		}
		return result;
	}

};

#endif
