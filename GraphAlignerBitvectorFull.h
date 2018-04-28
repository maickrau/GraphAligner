#ifndef GraphAlignerBitvectorFull_h
#define GraphAlignerBitvectorFull_h

#include "GraphAlignerCommon.h"
#include "GraphAlignerBitvectorCommon.h"
#include "WordSlice.h"
#include "NodeSlice.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorFull
{
private:
	using WordSlice = typename WordContainer<LengthType, ScoreType, Word>::Slice;
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using BV = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>;
	using EqVector = typename BV::EqVector;
	using Params = typename Common::Params;
	const Params params;
	const std::vector<std::vector<size_t>> componentOrder;
	const std::vector<size_t> belongsToComponent;
	class EdgeWithPriority
	{
	public:
		EdgeWithPriority(LengthType targetNode, ScoreType priority, bool full) :
		targetNode(targetNode),
		priority(priority),
		full(full)
		{
		}
		LengthType targetNode;
		ScoreType priority;
		bool full;
		bool operator>(const EdgeWithPriority& other) const
		{
			return priority > other.priority;
		}
	};
public:
	GraphAlignerBitvectorFull(const Params& params, const std::vector<std::vector<size_t>>& componentOrder, const std::vector<size_t>& belongsToComponent) :
	params(params),
	componentOrder(componentOrder),
	belongsToComponent(belongsToComponent)
	{
	}

	void verifyCorrectness(const std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::string& sequence, const EqVector& EqV, size_t j) const
	{
		for (size_t w = 0; w < currentSlice.size(); w++)
		{
			int cellbycellValue = previousSlice[w].scoreEnd+1;
			char graphChar = params.graph.NodeSequences(w);
			auto EqVMatch = EqV.getEq(graphChar);
			bool match = Common::characterMatch(graphChar, sequence[j]);
			size_t node = params.graph.IndexToNode(w);
			assert((EqVMatch & 1) == match);
			if (j != 0)
			{
				if (w == params.graph.NodeStart(node))
				{
					for (auto neighbor : params.graph.inNeighbors[node])
					{
						size_t u = params.graph.NodeEnd(neighbor)-1;
						cellbycellValue = std::min(cellbycellValue, previousSlice[u].scoreEnd+(match?0:1));
						cellbycellValue = std::min(cellbycellValue, currentSlice[u].getValue(0)+1);
					}
				}
				else
				{
					cellbycellValue = std::min(cellbycellValue, previousSlice[w-1].scoreEnd+(match?0:1));
					cellbycellValue = std::min(cellbycellValue, currentSlice[w-1].getValue(0)+1);
				}
				assert(currentSlice[w].getValue(0) == cellbycellValue);
			}
			else
			{
				assert(currentSlice[w].getValue(0) == (match ? 0 : 1));
			}
			for (size_t cell = 1; cell < WordConfiguration<Word>::WordSize && j+cell < sequence.size(); cell++)
			{
				match = Common::characterMatch(sequence[j+cell], graphChar);
				cellbycellValue = currentSlice[w].getValue(cell-1)+1;
				assert(((EqVMatch >> cell) & 1) == match);
				if (w == params.graph.NodeStart(node))
				{
					for (auto neighbor : params.graph.inNeighbors[node])
					{
						size_t u = params.graph.NodeEnd(neighbor)-1;
						cellbycellValue = std::min(cellbycellValue, currentSlice[u].getValue(cell-1)+(match?0:1));
						cellbycellValue = std::min(cellbycellValue, currentSlice[u].getValue(cell)+1);
					}
				}
				else
				{
					cellbycellValue = std::min(cellbycellValue, currentSlice[w-1].getValue(cell-1)+(match?0:1));
					cellbycellValue = std::min(cellbycellValue, currentSlice[w-1].getValue(cell)+1);
				}
				assert(currentSlice[w].getValue(cell) == cellbycellValue);
			}
		}
	}

	ScoreType alignAndGetScoreAcyclic(const std::string& sequence) const
	{
		std::vector<WordSlice> currentSlice;
		std::vector<WordSlice> previousSlice;
		currentSlice.resize(params.graph.NodeSequencesSize(), { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::WordSize, 0, true });
		previousSlice.resize(params.graph.NodeSequencesSize(), { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, 0, 0, true });
		int slices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		for (int slice = 0; slice < slices; slice++)
		{
			size_t j = slice * WordConfiguration<Word>::WordSize;
			EqVector EqV = BV::getEqVector(sequence, j);
			for (auto component : componentOrder)
			{
				assert(component.size() == 1);
				calculateNodeAcyclic(currentSlice, previousSlice, sequence, EqV, j, component[0]);
			}
			// verifyCorrectness(currentSlice, previousSlice, sequence, EqV, j);
			std::swap(previousSlice, currentSlice);
		}
		size_t extra = slices * WordConfiguration<Word>::WordSize - sequence.size();
		size_t scorepos = WordConfiguration<Word>::WordSize - extra;
		assert(scorepos < WordConfiguration<Word>::WordSize);

		ScoreType minScore = std::numeric_limits<ScoreType>::max();
		for (size_t i = 0; i < previousSlice.size(); i++)
		{
			minScore = std::min(minScore, previousSlice[i].getValue(scorepos));
		}
		return minScore;
	}

	ScoreType alignAndGetScore(const std::string& sequence) const
	{
		std::vector<WordSlice> currentSlice;
		std::vector<WordSlice> previousSlice;
		currentSlice.resize(params.graph.NodeSequencesSize(), { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::WordSize, 0, true });
		previousSlice.resize(params.graph.NodeSequencesSize(), { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, 0, 0, true });
		int slices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		for (int slice = 0; slice < slices; slice++)
		{
			size_t j = slice * WordConfiguration<Word>::WordSize;
			EqVector EqV = BV::getEqVector(sequence, j);
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
						calculateNodeAcyclic(currentSlice, previousSlice, sequence, EqV, j, component[0]);
						currentSlice[params.graph.NodeEnd(component[0])-1].calcMinScore();
						currentSlice[params.graph.NodeEnd(component[0])-1].sliceExists = true;
						continue;
					}
				}
				calculateCyclicComponent(sequence, component, EqV, currentSlice, previousSlice, j);
			}
			// verifyCorrectness(currentSlice, previousSlice, sequence, EqV, j);
			std::swap(previousSlice, currentSlice);
			currentSlice.assign(currentSlice.size(), { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::WordSize, 0, true });
		}
		size_t extra = slices * WordConfiguration<Word>::WordSize - sequence.size();
		size_t scorepos = WordConfiguration<Word>::WordSize - extra;
		assert(scorepos < WordConfiguration<Word>::WordSize);

		ScoreType minScore = std::numeric_limits<ScoreType>::max();
		for (size_t i = 0; i < previousSlice.size(); i++)
		{
			minScore = std::min(minScore, previousSlice[i].getValue(scorepos));
		}
		return minScore;
	}
private:

	void calculateCyclicComponent(const std::string& sequence, const std::vector<size_t>& nodes, const EqVector& EqV, std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, size_t j) const
	{
		std::priority_queue<EdgeWithPriority, std::vector<EdgeWithPriority>, std::greater<EdgeWithPriority>> calculables;
		size_t stillInQueues = 0;
		for (auto node : nodes)
		{
			auto start = params.graph.NodeStart(node);
			ScoreType startPriority = std::numeric_limits<ScoreType>::max();
			for (size_t i = 0; i < params.graph.NodeLength(node); i++)
			{
				currentSlice[start+i] = { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousSlice[start+i].scoreEnd+WordConfiguration<Word>::WordSize, previousSlice[start+i].scoreEnd, true };
				startPriority = std::min(startPriority, previousSlice[start+i].scoreEnd);
			}
			calculables.emplace(node, startPriority, true);
			for (auto neighbor : params.graph.inNeighbors[node])
			{
				if (belongsToComponent[neighbor] < belongsToComponent[node])
				{
					auto u = params.graph.NodeEnd(neighbor)-1;
					auto slice = currentSlice[u];
					assert(slice.sliceExists);
					assert(slice.minScore >= 0);
					assert(slice.minScore != std::numeric_limits<ScoreType>::max());
					auto graphChar = params.graph.NodeSequences(start);
					Word Eq = EqV.getEq(graphChar);
					WordSlice newSlice = BV::getNextSliceFullBand(Eq, currentSlice[u], (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[u]);
					if (newSlice.scoreBeforeStart > previousSlice[start].scoreEnd)
					{
						auto vertical = getSourceSliceFromScore(previousSlice[start].scoreEnd);
						newSlice = newSlice.mergeWithVertical(vertical);
						newSlice.scoreBeforeExists = true;
					}
					if (currentSlice[start].sliceExists)
					{
						newSlice = newSlice.mergeWith(currentSlice[start]);
					}
					newSlice.sliceExists = true;
					currentSlice[start] = newSlice;
				}
			}
		}
		while (calculables.size() > 0)
		{
			auto top = calculables.top();
			calculables.pop();
			stillInQueues--;
			auto start = params.graph.NodeStart(top.targetNode);
			auto end = params.graph.NodeEnd(top.targetNode);
			bool calculatedToEnd = true;
			for (size_t w = start+1; w < end; w++)
			{
				auto graphChar = params.graph.NodeSequences(w);
				Word Eq = EqV.getEq(graphChar);
				WordSlice oldSlice = currentSlice[w];
				WordSlice newSlice = BV::getNextSliceFullBand(Eq, currentSlice[w-1], (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[w-1]);
				if (newSlice.scoreBeforeStart > previousSlice[w].scoreEnd)
				{
					auto vertical = getSourceSliceFromScore(previousSlice[w].scoreEnd);
					newSlice = newSlice.mergeWithVertical(vertical);
					newSlice.scoreBeforeExists = true;
				}
				if (oldSlice.sliceExists)
				{
					BV::assertSliceCorrectness(oldSlice, previousSlice[w], true);
					newSlice = newSlice.mergeWith(oldSlice);
					if (!top.full)
					{
						auto newPriority = newSlice.changedMinScore(oldSlice);
						if (newPriority == std::numeric_limits<ScoreType>::max())
						{
							calculatedToEnd = false;
							break;
						}
					}
				}
				newSlice.sliceExists = true;
				currentSlice[w] = newSlice;
			}
			if (calculatedToEnd)
			{
				for (auto neighbor : params.graph.outNeighbors[top.targetNode])
				{
					if (belongsToComponent[neighbor] == belongsToComponent[top.targetNode])
					{
						auto u = params.graph.NodeStart(neighbor);
						auto graphChar = params.graph.NodeSequences(u);
						Word Eq = EqV.getEq(graphChar);
						WordSlice oldSlice = currentSlice[u];
						WordSlice newSlice = BV::getNextSliceFullBand(Eq, currentSlice[end-1], (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[end-1]);
						if (newSlice.scoreBeforeStart > previousSlice[u].scoreEnd)
						{
							auto vertical = getSourceSliceFromScore(previousSlice[u].scoreEnd);
							newSlice = newSlice.mergeWithVertical(vertical);
							newSlice.scoreBeforeExists = true;
						}
						ScoreType newPriority;
						BV::assertSliceCorrectness(oldSlice, previousSlice[u], true);
						newSlice = newSlice.mergeWith(oldSlice);
						newSlice.calcMinScore();
						newPriority = newSlice.changedMinScore(oldSlice);
						BV::assertSliceCorrectness(newSlice, previousSlice[u], true);
						if (newPriority == std::numeric_limits<ScoreType>::max())
						{
							calculatedToEnd = false;
							continue;
						}
						newSlice.sliceExists = true;
						currentSlice[u] = newSlice;
						calculables.emplace(neighbor, newPriority, false);
						stillInQueues++;
					}
				}
			}
		}
	}

	WordSlice getSourceSliceFromScore(ScoreType previousScore) const
	{
		WordSlice result { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousScore+WordConfiguration<Word>::WordSize, previousScore, false };
		result.sliceExists = true;
		return result;
	}

	WordSlice getStartSliceAcyclic(const std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::string& sequence, const EqVector& EqV, size_t j, size_t node, size_t start, const bool previousEq) const
	{
		if (params.graph.inNeighbors[node].size() == 0)
		{
			if (j == 0)
			{
				auto result = getSourceSliceFromScore(0);
				if (sequence[0] == params.graph.NodeSequences(params.graph.NodeStart(node)))
				{
					assert(result.VP & 1);
					result.VP &= ~1;
					result.scoreEnd -= 1;
				}
				return result;
			}
			else
			{
				return getSourceSliceFromScore(previousSlice[start].scoreEnd);
			}
		}
		char graphChar = params.graph.NodeSequences(start);
		Word Eq = EqV.getEq(graphChar);
		if (params.graph.inNeighbors[node].size() == 1)
		{
			size_t u = params.graph.NodeEnd(params.graph.inNeighbors[node][0])-1;
			return BV::getNextSliceFullBand(Eq, currentSlice[u], (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[u]);
		}
		WordSlice result;
		WordSlice up;
		up = previousSlice[start];
		bool foundOne = false;
		for (auto neighbor : params.graph.inNeighbors[node])
		{
			WordSlice previous;
			WordSlice previousUp;
			size_t u = params.graph.NodeEnd(neighbor)-1;
			previous = currentSlice[u];
			previousUp = previousSlice[u];
			assert(previous.scoreBeforeStart == previousUp.scoreEnd);
			BV::assertSliceCorrectness(previous, previousUp, true);
			auto resultHere = BV::getNextSliceFullBand(Eq, previous, previousEq, previousUp);
			if (!foundOne)
			{
				result = resultHere;
				foundOne = true;
			}
			else
			{
				result = result.mergeWith(resultHere);
			}
		}
		assert(foundOne);
		return result;
	}

	void calculateNodeAcyclic(std::vector<WordSlice>& currentSlice, const std::vector<WordSlice>& previousSlice, const std::string& sequence, const EqVector& EqV, size_t j, size_t node) const
	{
		auto start = params.graph.NodeStart(node);
		auto end = params.graph.NodeEnd(node);
		auto length = end-start;
		currentSlice[start] = getStartSliceAcyclic(currentSlice, previousSlice, sequence, EqV, j, node, start, (j == 0) || (j > 0 && params.graph.NodeSequences(start) == sequence[j-1]));
		currentSlice[start].scoreBeforeExists = true;
		if (currentSlice[start].scoreBeforeStart > previousSlice[start].scoreEnd)
		{
			auto vertical = getSourceSliceFromScore(previousSlice[start].scoreEnd);
			currentSlice[start] = currentSlice[start].mergeWithVertical(vertical);
			currentSlice[start].scoreBeforeExists = true;
		}
		BV::assertSliceCorrectness(currentSlice[start], previousSlice[start], true);
		for (size_t w = start+1; w < end; w++)
		{
			char graphChar = params.graph.NodeSequences(w);
			Word Eq = EqV.getEq(graphChar);
			currentSlice[w] = BV::getNextSliceFullBand(Eq, currentSlice[w-1], (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[w-1]);
			if (currentSlice[w].scoreBeforeStart > previousSlice[w].scoreEnd)
			{
				auto vertical = getSourceSliceFromScore(previousSlice[w].scoreEnd);
				currentSlice[w] = currentSlice[w].mergeWithVertical(vertical);
			}
			currentSlice[w].scoreBeforeExists = true;
			BV::assertSliceCorrectness(currentSlice[w], previousSlice[w], true);
		}
	}

};

#endif
