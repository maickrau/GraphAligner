#ifndef GraphAlignerCellbycellFull_h
#define GraphAlignerCellbycellFull_h

#include "GraphAlignerCommon.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerCellbycellFull
{
private:
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using Params = typename Common::Params;
	const Params params;
	const std::vector<std::vector<size_t>> componentOrder;
public:
	GraphAlignerCellbycellFull(const Params& params, const std::vector<std::vector<size_t>>& componentOrder) :
	params(params),
	componentOrder(componentOrder)
	{
	}

	ScoreType alignAndGetScoreAcyclic(const std::string& sequence) const
	{
		std::vector<ScoreType> currentSlice;
		std::vector<ScoreType> previousSlice;
		currentSlice.resize(params.graph.NodeSequencesSize(), 0);
		previousSlice.resize(params.graph.NodeSequencesSize(), 0);
		for (size_t w = 0; w < previousSlice.size(); w++)
		{
			char graphChar = params.graph.NodeSequences(w);
			bool match = Common::characterMatch(sequence[0], graphChar);
			previousSlice[w] = match ? 0 : 1;
		}
		for (size_t j = 1; j < sequence.size(); j++)
		{
			for (auto component : componentOrder)
			{
				assert(component.size() == 1);
				calculateNodeAcyclic(currentSlice, previousSlice, sequence, j, component[0]);
			}
			std::swap(currentSlice, previousSlice);
		}

		ScoreType minScore = std::numeric_limits<ScoreType>::max();
		for (size_t i = 0; i < previousSlice.size(); i++)
		{
			minScore = std::min(minScore, previousSlice[i]);
		}
		return minScore;
	}

private:

	void calculateNodeAcyclic(std::vector<ScoreType>& currentSlice, const std::vector<ScoreType>& previousSlice, const std::string& sequence, size_t j, size_t node) const
	{
		auto start = params.graph.NodeStart(node);
		auto end = params.graph.NodeEnd(node);
		currentSlice[start] = previousSlice[start]+1;
		char graphChar = params.graph.NodeSequences(start);
		bool match = Common::characterMatch(sequence[j], graphChar);
		for (auto neighbor : params.graph.inNeighbors[node])
		{
			size_t u = params.graph.NodeEnd(neighbor)-1;
			currentSlice[start] = std::min(currentSlice[start], std::min(currentSlice[u]+1, previousSlice[u]+(match?0:1)));
		}
		for (size_t w = start+1; w < end; w++)
		{
			graphChar = params.graph.NodeSequences(w);
			match = Common::characterMatch(sequence[j], graphChar);
			currentSlice[w] = std::min(currentSlice[w-1]+1, std::min(previousSlice[w]+1, previousSlice[w-1]+(match?0:1)));
		}
	}

};

#endif
