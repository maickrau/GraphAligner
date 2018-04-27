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
public:
	GraphAlignerBitvectorFull(const Params& params, const std::vector<std::vector<size_t>>& componentOrder) :
	params(params),
	componentOrder(componentOrder)
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
			std::swap(currentSlice, previousSlice);
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

	// ScoreType alignAndGetScore(const std::string& sequence) const
	// {
	// 	std::vector<WordSlice> currentSlice;
	// 	std::vector<WordSlice> previousSlice;
	// 	currentSlice.resize(params.graph.NodeSequencesSize() { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::WordSize, 0, true });
	// 	previousSlice.resize(params.graph.NodeSequencesSize() { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, 0, 0, true });
	// 	ArrayPriorityQueue<NodeWithPriority> calculableQueue { sequence.size() };
	// 	int slices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
	// 	for (int slice = 0; slice < slices.size(); slice++)
	// 	{
	// 		size_t j = slice * WordConfiguration<Word>::WordSize;
	// 		EqVector EqV = BV::getEqVector(sequence, j);

	// 		while (calculableQueue.size() > 0)
	// 		{
				
	// 		}
	// 		std::swap(previousSlice, currentSlice);
	// 	}
	// }
private:

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
			return BV::getNextSlice(Eq, currentSlice[u], true, true, true, (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[u], previousSlice[start]);
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
			auto resultHere = BV::getNextSlice(Eq, previous, true, true, true, previousEq, previousUp, currentSlice[start]);
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
		if (currentSlice[start].scoreBeforeStart > previousSlice[start].scoreEnd)
		{
			auto vertical = getSourceSliceFromScore(previousSlice[start].scoreEnd);
			currentSlice[start] = currentSlice[start].mergeWithVertical(vertical);
			currentSlice[start].scoreBeforeExists = true;
		}
		for (size_t w = start+1; w < end; w++)
		{
			char graphChar = params.graph.NodeSequences(w);
			Word Eq = EqV.getEq(graphChar);
			currentSlice[w] = BV::getNextSlice(Eq, currentSlice[w-1], true, true, true, (j == 0) || (j > 0 && graphChar == sequence[j-1]), previousSlice[w-1], previousSlice[w]);
			if (currentSlice[w].scoreBeforeStart > previousSlice[w].scoreEnd)
			{
				auto vertical = getSourceSliceFromScore(previousSlice[w].scoreEnd);
				currentSlice[w] = currentSlice[w].mergeWithVertical(vertical);
			}
			currentSlice[w].scoreBeforeExists = true;
		}
	}

};

#endif
