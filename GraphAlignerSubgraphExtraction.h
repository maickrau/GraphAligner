#ifndef GraphAlignerSubgraphExtraction_h
#define GraphAlignerSubgraphExtraction_h

#include <unordered_map>
#include "AlignmentGraph.h"

template <typename LengthType, typename ScoreType, typename Word>
class SubgraphExtractor
{
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	using SeedHit = typename Common::SeedHit;
public:
	static AlignmentGraph ExtractSubgraph(AlignerGraphsizedState& state, const AlignmentGraph& original, const std::vector<SeedHit>& seedHits, const size_t maxDist)
	{
		return extract(state, original, seedHits, maxDist);
	}
private:
	static AlignmentGraph extract(AlignerGraphsizedState& state, const AlignmentGraph& original, const std::vector<SeedHit>& seedHits, const size_t maxDist)
	{
		std::vector<std::vector<size_t>> propagateFw;
		std::vector<std::vector<size_t>> propagateBw;
		propagateFw.resize(maxDist);
		propagateBw.resize(maxDist);
		std::vector<size_t> check;
		for (auto seedHit : seedHits)
		{
			LengthType forwardNodeId;
			if (seedHit.reverse)
			{
				forwardNodeId = seedHit.nodeID * 2 + 1;
			}
			else
			{
				forwardNodeId = seedHit.nodeID * 2;
			}
			auto revPos = original.GetReversePosition(forwardNodeId, seedHit.nodeOffset);
			auto nodeIndex = original.GetUnitigNode(forwardNodeId, seedHit.nodeOffset);
			auto revNodeIndex = original.GetUnitigNode(revPos.first, revPos.second);
			assert(nodeIndex < state.forwardNodeDist.size());
			assert(state.forwardNodeDist[nodeIndex] == std::numeric_limits<size_t>::max() || state.forwardNodeDist[nodeIndex] == 0);
			assert(state.backwardNodeDist[nodeIndex] == std::numeric_limits<size_t>::max() || state.backwardNodeDist[nodeIndex] == 0);
			assert(revNodeIndex < state.forwardNodeDist.size());
			assert(state.forwardNodeDist[revNodeIndex] == std::numeric_limits<size_t>::max() || state.forwardNodeDist[revNodeIndex] == 0);
			assert(state.backwardNodeDist[revNodeIndex] == std::numeric_limits<size_t>::max() || state.backwardNodeDist[revNodeIndex] == 0);
			propagateFw[0].push_back(nodeIndex);
			propagateBw[0].push_back(nodeIndex);
			propagateBw[0].push_back(revNodeIndex);
			propagateFw[0].push_back(revNodeIndex);
		}

		for (size_t dist = 0; dist < maxDist; dist++)
		{
			for (auto fw : propagateFw[dist])
			{
				if (state.forwardNodeDist[fw] <= dist) continue;
				state.forwardNodeDist[fw] = dist;
				check.push_back(fw);
				if (dist + original.NodeLength(fw) >= maxDist) continue;
				size_t newDist = dist + original.NodeLength(fw);
				assert(newDist > dist);
				for (auto neighbor : original.outNeighbors[fw])
				{
					propagateFw[newDist].push_back(neighbor);
				}
			}
			for (auto bw : propagateBw[dist])
			{
				if (state.backwardNodeDist[bw] <= dist) continue;
				state.backwardNodeDist[bw] = dist;
				check.push_back(bw);
				if (dist + original.NodeLength(bw) >= maxDist) continue;
				size_t newDist = dist + original.NodeLength(bw);
				assert(newDist > dist);
				for (auto neighbor : original.inNeighbors[bw])
				{
					propagateBw[newDist].push_back(neighbor);
				}
			}
		}

		size_t subgraphSize = 0;
		std::unordered_map<size_t, size_t> nodeMapping;
		for (auto node : check)
		{
			if (nodeMapping.count(node) == 1) continue;
			if (state.forwardNodeDist[node] != std::numeric_limits<size_t>::max() && state.backwardNodeDist[node] != std::numeric_limits<size_t>::max() && state.forwardNodeDist[node] + state.backwardNodeDist[node] <= maxDist)
			{
				nodeMapping[node] = subgraphSize;
				subgraphSize++;
			}
			state.forwardNodeDist[node] = std::numeric_limits<size_t>::max();
			state.backwardNodeDist[node] = std::numeric_limits<size_t>::max();
		}

		return original.GetSubgraph(nodeMapping);
	}
};

#endif
