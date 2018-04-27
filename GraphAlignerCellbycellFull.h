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
	const std::vector<size_t> belongsToComponent;
public:
	GraphAlignerCellbycellFull(const Params& params, const std::vector<std::vector<size_t>>& componentOrder, const std::vector<size_t>& belongsToComponent) :
	params(params),
	componentOrder(componentOrder),
	belongsToComponent(belongsToComponent)
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
			// verifyCorrectness(currentSlice, previousSlice, sequence[j]);
			std::swap(currentSlice, previousSlice);
		}

		ScoreType minScore = std::numeric_limits<ScoreType>::max();
		for (size_t i = 0; i < previousSlice.size(); i++)
		{
			minScore = std::min(minScore, previousSlice[i]);
		}
		return minScore;
	}

	void verifyCorrectness(const std::vector<ScoreType>& currentSlice, const std::vector<ScoreType>& previousSlice, const char sequenceChar) const
	{
		for (size_t w = 0; w < currentSlice.size(); w++)
		{
			auto graphChar = params.graph.NodeSequences(w);
			bool match = Common::characterMatch(sequenceChar, graphChar);
			auto node = params.graph.IndexToNode(w);
			ScoreType foundScore = previousSlice[w]+1;
			if (w == params.graph.NodeStart(node))
			{
				for (auto neighbor : params.graph.inNeighbors[node])
				{
					auto u = params.graph.NodeEnd(neighbor)-1;
					foundScore = std::min(foundScore, currentSlice[u]+1);
					foundScore = std::min(foundScore, previousSlice[u]+(match?0:1));
				}
			}
			else
			{
				auto u = w-1;
				foundScore = std::min(foundScore, currentSlice[u]+1);
				foundScore = std::min(foundScore, previousSlice[u]+(match?0:1));
			}
			assert(currentSlice[w] == foundScore);
		}
	}

	ScoreType alignAndGetScore(const std::string& sequence) const
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
			currentSlice.assign(currentSlice.size(), sequence.size());
			for (auto component : componentOrder)
			{
				if (component.size() == 1)
				{
					bool isAcyclic = true;
					for (auto neighbor : params.graph.inNeighbors[component[0]])
					{
						if (neighbor == component[0])
						{
							isAcyclic = false;
							break;
						}
					}
					if (isAcyclic)
					{
						calculateNodeAcyclic(currentSlice, previousSlice, sequence, j, component[0]);
						continue;
					}
				}
				calculateCyclicComponent(sequence, component, currentSlice, previousSlice, j);
			}
			// verifyCorrectness(currentSlice, previousSlice, sequence[j]);
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

	void calculateCyclicComponent(const std::string& sequence, const std::vector<size_t>& nodes, std::vector<ScoreType>& currentSlice, const std::vector<ScoreType>& previousSlice, size_t j) const
	{
		for (auto node : nodes)
		{
			calculateNodeAcyclic(currentSlice, previousSlice, sequence, j, node);
		}
		for (auto node : nodes)
		{
			auto start = params.graph.NodeStart(node);
			for (auto neighbor : params.graph.inNeighbors[node])
			{
				auto u = params.graph.NodeEnd(neighbor)-1;
				if (currentSlice[start] > currentSlice[u] + 1)
				{
					recurseHorizontalScores(currentSlice, node, start, currentSlice[u]+1);
				}
			}
			for (LengthType w = start+1; w < params.graph.NodeEnd(node); w++)
			{
				if (currentSlice[w] > currentSlice[w-1]+1)
				{
					recurseHorizontalScores(currentSlice, node, w, currentSlice[w-1]+1);
				}
			}
		}
	}

	void recurseHorizontalScores(std::vector<ScoreType>& currentSlice, size_t node, size_t start, ScoreType newScore) const
	{
		assert(currentSlice[start] > newScore);
		currentSlice[start] = newScore;
		newScore++;
		for (LengthType w = start+1; w < params.graph.NodeEnd(node); w++)
		{
			if (currentSlice[w] <= newScore) return;
			currentSlice[w] = newScore;
			newScore++;
		}
		for (auto neighbor : params.graph.outNeighbors[node])
		{
			if (belongsToComponent[neighbor] == belongsToComponent[node])
			{
				auto u = params.graph.NodeStart(neighbor);
				if (currentSlice[u] > newScore)
				{
					recurseHorizontalScores(currentSlice, neighbor, u, newScore);
				}
			}
		}
	}

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
