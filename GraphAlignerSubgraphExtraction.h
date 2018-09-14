#ifndef GraphAlignerSubgraphExtraction_h
#define GraphAlignerSubgraphExtraction_h

#include <unordered_map>
#include "GraphAlignerWrapper.h"
#include "AlignmentGraph.h"

template <typename LengthType, typename ScoreType, typename Word>
class SubgraphExtractor
{
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
public:
	static void ExtractSubgraph(AlignerGraphsizedState& state, const AlignmentGraph& original, const std::vector<SeedHit>& seedHits, size_t readSize)
	{
		extract(state, original, seedHits, readSize);
	}
private:
	static void extract(AlignerGraphsizedState& state, const AlignmentGraph& original, std::vector<SeedHit> seedHits, size_t readSize)
	{
		std::sort(seedHits.begin(), seedHits.end(), [](const SeedHit& left, const SeedHit& right){ return left.seqPos < right.seqPos; });
		std::vector<std::vector<size_t>> propagateFw;
		std::vector<std::vector<size_t>> propagateBw;
		// size_t maxLen = 0;
		// if (seedHits.size() < 5)
		// {
		// 	maxLen = seedHits.back().seqPos - seedHits[0].seqPos;
		// }
		// else
		// {
		// 	for (size_t i = 0; i < seedHits.size() - 5; i++)
		// 	{
		// 		maxLen = std::max(maxLen, seedHits[i+4].seqPos - seedHits[i].seqPos);
		// 	}
		// }
		size_t maxLen = readSize / seedHits.size();
		// size_t maxLen = readSize;
		if (maxLen == 0) maxLen = 1;
		addNodes(state, original, seedHits, propagateBw, propagateFw, maxLen);
		// if (seedHits.size() < 5)
		// {
		// 	addNodes(state, original, seedHits, propagateBw, propagateFw, seedHits.back().seqPos - seedHits[0].seqPos);
		// }
		// else
		// {
		// 	for (size_t i = 0; i <= seedHits.size() - 5; i++)
		// 	{
		// 		addNodes(state, original, {seedHits.begin()+i, seedHits.begin()+i+5}, propagateBw, propagateFw, seedHits[i+4].seqPos - seedHits[i].seqPos);
		// 	}
		// }
	}

	static void addNodes(AlignerGraphsizedState& state, const AlignmentGraph& original, std::vector<SeedHit> seedHits, std::vector<std::vector<size_t>> propagateBw, std::vector<std::vector<size_t>> propagateFw, const size_t maxDist)
	{
		assert(maxDist > 0);
		if (propagateFw.size() < maxDist) propagateFw.resize(maxDist);
		if (propagateBw.size() < maxDist) propagateBw.resize(maxDist);
		std::vector<size_t> check;
		for (auto seedHit : seedHits)
		{
			LengthType forwardNodeId;
			LengthType backwardNodeId;
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
			auto pos = original.unitigReidLookup.at(forwardNodeId);
			auto revPos = original.unitigReidLookup.at(backwardNodeId);
			assert(pos < state.forwardReidDist.size());
			assert(pos < state.backwardReidDist.size());
			assert(state.forwardReidDist[pos] == std::numeric_limits<size_t>::max() || state.forwardReidDist[pos] == 0);
			assert(state.backwardReidDist[pos] == std::numeric_limits<size_t>::max() || state.backwardReidDist[pos] == 0);
			assert(revPos < state.forwardReidDist.size());
			assert(revPos < state.backwardReidDist.size());
			assert(state.forwardReidDist[revPos] == std::numeric_limits<size_t>::max() || state.forwardReidDist[revPos] == 0);
			assert(state.backwardReidDist[revPos] == std::numeric_limits<size_t>::max() || state.backwardReidDist[revPos] == 0);
			state.forwardReidDist[pos] = 0;
			state.backwardReidDist[pos] = 0;
			state.forwardReidDist[revPos] = 0;
			state.backwardReidDist[revPos] = 0;
			if (!state.subgraph[pos])
			{
				state.subgraph[pos] = true;
				state.subgraphReids.push_back(pos);
			}
			if (!state.subgraph[revPos])
			{
				state.subgraph[revPos] = true;
				state.subgraphReids.push_back(revPos);
			}
			for (auto neighbor : original.unitigReidOutNeighbors[pos])
			{
				propagateFw[0].push_back(neighbor);
			}
			for (auto neighbor : original.unitigReidOutNeighbors[revPos])
			{
				propagateFw[0].push_back(neighbor);
			}
			for (auto neighbor : original.unitigReidInNeighbors[pos])
			{
				propagateBw[0].push_back(neighbor);
			}
			for (auto neighbor : original.unitigReidInNeighbors[revPos])
			{
				propagateBw[0].push_back(neighbor);
			}
		}

		for (size_t dist = 0; dist < maxDist; dist++)
		{
			for (auto fw : propagateFw[dist])
			{
				assert(fw < state.forwardReidDist.size());
				if (state.forwardReidDist[fw] <= dist) continue;
				state.forwardReidDist[fw] = dist;
				check.push_back(fw);
				size_t newDist = dist + original.ReidLength(fw) - original.DBGOverlap;
				if (newDist >= maxDist) continue;
				assert(newDist > dist);
				assert(newDist < propagateFw.size());
				for (auto neighbor : original.unitigReidOutNeighbors[fw])
				{
					propagateFw[newDist].push_back(neighbor);
				}
			}
			for (auto bw : propagateBw[dist])
			{
				assert(bw < state.backwardReidDist.size());
				if (state.backwardReidDist[bw] <= dist) continue;
				state.backwardReidDist[bw] = dist;
				check.push_back(bw);
				size_t newDist = dist + original.ReidLength(bw) - original.DBGOverlap;
				if (newDist >= maxDist) continue;
				assert(newDist > dist);
				assert(newDist < propagateBw.size());
				for (auto neighbor : original.unitigReidInNeighbors[bw])
				{
					propagateBw[newDist].push_back(neighbor);
				}
			}
			propagateFw[dist].clear();
			propagateBw[dist].clear();
		}

		for (auto node : check)
		{
			assert(node < state.subgraph.size());
			if (state.subgraph[node])
			{
				state.forwardReidDist[node] = std::numeric_limits<size_t>::max();
				state.backwardReidDist[node] = std::numeric_limits<size_t>::max();
				continue;
			}
			if (state.forwardReidDist[node] != std::numeric_limits<size_t>::max() && state.backwardReidDist[node] != std::numeric_limits<size_t>::max() && state.forwardReidDist[node] + state.backwardReidDist[node] <= 2 * maxDist)
			{
				state.subgraph[node] = true;
				state.subgraphReids.push_back(node);
			}
			state.forwardReidDist[node] = std::numeric_limits<size_t>::max();
			state.backwardReidDist[node] = std::numeric_limits<size_t>::max();
		}

#ifndef NDEBUG
		for (auto seed : seedHits)
		{
			LengthType forwardNodeId;
			LengthType backwardNodeId;
			if (seed.reverse)
			{
				forwardNodeId = seed.nodeID * 2 + 1;
				backwardNodeId = seed.nodeID * 2;
			}
			else
			{
				forwardNodeId = seed.nodeID * 2;
				backwardNodeId = seed.nodeID * 2 + 1;
			}
			assert(state.subgraph[forwardNodeId]);
			assert(state.subgraph[backwardNodeId]);
		}
#endif

	}
};

#endif
