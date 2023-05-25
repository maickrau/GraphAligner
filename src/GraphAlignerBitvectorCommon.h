#ifndef GraphAlignerBitvectorCommon_h
#define GraphAlignerBitvectorCommon_h

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "AlignmentGraph.h"
#include "NodeSlice.h"
#include "CommonUtils.h"
#include "GraphAlignerWrapper.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"
#include "GraphAlignerCommon.h"
#include "ArrayPriorityQueue.h"

#ifndef NDEBUG
thread_local int debugLastRowMinScore;
#endif

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorCommon
{
private:
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using AlignerGraphsizedState = typename Common::AlignerGraphsizedState;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
	using Trace = typename Common::Trace;
	using OnewayTrace = typename Common::OnewayTrace;
	using EdgeWithPriority = typename Common::EdgeWithPriority;
public:
	using WordSlice = decltype(NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem::startSlice);

	class EqVector
	{
	public:
		EqVector(Word BA, Word BT, Word BC, Word BG)
		{
			masks[0] = BA;
			masks[1] = BC;
			masks[2] = BG;
			masks[3] = BT;
		}
		Word getEqI(AlignmentGraph::AmbiguousChunkSequence eq) const
		{
			assert((eq.A | eq.C | eq.G | eq.T) & 1);
			Word result = 0;
			if (eq.A & 1) result |= A();
			if (eq.C & 1) result |= C();
			if (eq.G & 1) result |= G();
			if (eq.T & 1) result |= T();
			return result;
		}
		Word getEqC(char c) const
		{
			switch(c)
			{
				case 'A':
				case 'a':
					return A();
				case 'C':
				case 'c':
					return C();
				case 'G':
				case 'g':
					return G();
				case 'T':
				case 't':
				case 'U':
				case 'u':
					return T();
				case 'R':
				case 'r':
					return A() | G();
				case 'Y':
				case 'y':
					return C() | T();
				case 'S':
				case 's':
					return G() | C();
				case 'W':
				case 'w':
					return A() | T();
				case 'K':
				case 'k':
					return G() | T();
				case 'M':
				case 'm':
					return A() | C();
				case 'B':
				case 'b':
					return C() | G() | T();
				case 'D':
				case 'd':
					return A() | G() | T();
				case 'H':
				case 'h':
					return A() | C() | T();
				case 'V':
				case 'v':
					return A() | C() | G();
				case 'N':
				case 'n':
					return A() | C() | G() | T();
				default:
					assert(false);
			}
			assert(false);
			return 0;
		}
		Word getEqI(size_t i) const
		{
			assert(i < 4);
			return masks[i];
		}
		Word masks[4];
	private:
		Word A() const
		{
			return masks[0];
		}
		Word C() const
		{
			return masks[1];
		}
		Word G() const
		{
			return masks[2];
		}
		Word T() const
		{
			return masks[3];
		}
	};
	class DPSlice
	{
	public:
		DPSlice() :
		minScore(std::numeric_limits<ScoreType>::max()),
		minScoreNode(std::numeric_limits<LengthType>::max()),
		minScoreNodeOffset(std::numeric_limits<LengthType>::max()),
		maxExactEndposScore(std::numeric_limits<ScoreType>::min()),
		maxExactEndposNode(std::numeric_limits<LengthType>::max()),
		scoresVectorMap(),
		scores(),
		j(std::numeric_limits<LengthType>::max()),
		cellsProcessed(0),
		bandwidth(0),
		scoresNotValid(false),
		seedstartNodes()
#ifdef SLICEVERBOSE
		,nodesProcessed(0)
		,numCells(0)
#endif
		{}
		DPSlice(std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem>* vectorMap) :
		minScore(std::numeric_limits<ScoreType>::max()),
		minScoreNode(std::numeric_limits<LengthType>::max()),
		minScoreNodeOffset(std::numeric_limits<LengthType>::max()),
		maxExactEndposScore(std::numeric_limits<ScoreType>::min()),
		maxExactEndposNode(std::numeric_limits<LengthType>::max()),
		scoresVectorMap(vectorMap),
		scores(),
		j(std::numeric_limits<LengthType>::max()),
		cellsProcessed(0),
		bandwidth(0),
		scoresNotValid(false)
#ifdef SLICEVERBOSE
		,nodesProcessed(0)
		,numCells(0)
#endif
		{}
		ScoreType minScore;
		LengthType minScoreNode;
		LengthType minScoreNodeOffset;
		ScoreType maxExactEndposScore;
		LengthType maxExactEndposNode;
		NodeSlice<LengthType, ScoreType, Word, true> scoresVectorMap;
		NodeSlice<LengthType, ScoreType, Word, false> scores;
		LengthType j;
		size_t cellsProcessed;
		size_t bandwidth;
		bool scoresNotValid;
		std::unordered_set<size_t> seedstartNodes;
		phmap::flat_hash_map<size_t, ScoreType> nodeMaxExactEndposScore;
#ifdef SLICEVERBOSE
		size_t nodesProcessed;
		size_t numCells;
#endif
		DPSlice getMapSlice() const
		{
			DPSlice result;
			result.minScore = minScore;
			result.minScoreNode = minScoreNode;
			result.minScoreNodeOffset = minScoreNodeOffset;
			result.maxExactEndposNode = maxExactEndposNode;
			result.maxExactEndposScore = maxExactEndposScore;
			result.scores = scores;
			result.j = j;
			result.cellsProcessed = cellsProcessed;
			result.bandwidth = bandwidth;
			result.scoresNotValid = scoresNotValid;
			result.seedstartNodes = seedstartNodes;
			result.nodeMaxExactEndposScore = nodeMaxExactEndposScore;
#ifdef SLICEVERBOSE
			result.nodesProcessed = nodesProcessed;
			result.numCells = numCells;
#endif
			return result;
		}
	};
	class DPTable
	{
	public:
		DPTable() :
		slices()
		{}
		std::vector<DPSlice> slices;
	};
	class NodeCalculationResult
	{
	public:
		ScoreType minScore;
		LengthType minScoreNode;
		LengthType minScoreNodeOffset;
		LengthType maxExactEndposNode;
		ScoreType maxExactEndposScore;
		size_t cellsProcessed;
#ifdef SLICEVERBOSE
		size_t nodesProcessed;
#endif
	};


	GraphAlignerBitvectorCommon() = delete;

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static inline std::tuple<WordSlice, Word, Word> getNextSlice(Word Eq, WordSlice slice, Word hinP, Word hinN)
	{
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408

		Word Xv = Eq | slice.VN; //line 7
		Eq |= hinN; //between lines 7-8
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq; //line 8
		Word Ph = slice.VN | ~(Xh | slice.VP); //line 9
		Word Mh = slice.VP & Xh; //line 10
		Word tempMh = (Mh << 1) | hinN; //line 16 + between lines 16-17
		hinN = Mh >> (WordConfiguration<Word>::WordSize-1); //line 11
		Word tempPh = (Ph << 1) | hinP; //line 15 + between lines 16-17
		slice.VP = tempMh | ~(Xv | tempPh); //line 17
		hinP = Ph >> (WordConfiguration<Word>::WordSize-1); //line 13
		slice.VN = tempPh & Xv; //line 18
		slice.scoreEnd -= hinN; //line 12
		slice.scoreEnd += hinP; //line 14

		return std::make_tuple(slice, Ph, Mh);
	}

	static WordSlice flattenWordSlice(WordSlice slice, size_t row)
	{
		Word mask = ~(WordConfiguration<Word>::AllOnes << row);
		slice.scoreEnd -= WordConfiguration<Word>::popcount(slice.VP & ~mask);
		slice.scoreEnd += WordConfiguration<Word>::popcount(slice.VN & ~mask);
		slice.VP &= mask;
		slice.VN &= mask;
		return slice;
	}

	static EqVector getEqVector(const std::string& sequence, size_t j)
	{
		return getEqVector(std::string_view { sequence.data(), sequence.size() }, j);
	}

	static EqVector getEqVector(const std::string_view& sequence, size_t j)
	{
		Word BA = WordConfiguration<Word>::AllZeros;
		Word BT = WordConfiguration<Word>::AllZeros;
		Word BC = WordConfiguration<Word>::AllZeros;
		Word BG = WordConfiguration<Word>::AllZeros;
		for (int i = 0; i < WordConfiguration<Word>::WordSize && j+i < sequence.size(); i++)
		{
			Word mask = ((Word)1) << i;
			switch(sequence[j+i])
			{
				case 'a':
				case 'A':
					BA |= mask;
					break;
				case 'c':
				case 'C':
					BC |= mask;
					break;
				case 'g':
				case 'G':
					BG |= mask;
					break;
				case 't':
				case 'T':
					BT |= mask;
					break;
				default:
					if (Common::characterMatch(sequence[j+i], 'A')) BA |= mask;
					if (Common::characterMatch(sequence[j+i], 'C')) BC |= mask;
					if (Common::characterMatch(sequence[j+i], 'T')) BT |= mask;
					if (Common::characterMatch(sequence[j+i], 'G')) BG |= mask;
					break;
			}
		}
		assert((j + WordConfiguration<Word>::WordSize > sequence.size()) || (BA | BC | BT | BG) == WordConfiguration<Word>::AllOnes);
		assert((j + WordConfiguration<Word>::WordSize <= sequence.size()) || ((BA | BC | BT | BG) & (WordConfiguration<Word>::AllOnes << (WordConfiguration<Word>::WordSize - j + sequence.size()))) == WordConfiguration<Word>::AllZeros);
		EqVector EqV {BA, BT, BC, BG};
		return EqV;
	}

	static WordSlice getSeedSlice(size_t j, size_t seqLen, const Params& params)
	{
		assert(j % WordConfiguration<Word>::WordSize == 0);
		WordSlice result { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, 0 };
		ScoreType oldSeedScore = 0;
		if (j > 0)
		{
			oldSeedScore = ((ScoreType)(j-1+1)*100 + params.XscoreErrorCost / 2) / params.XscoreErrorCost;
		}
		for (size_t offset = 0; offset < WordConfiguration<Word>::WordSize; offset++)
		{
			if (j+offset >= seqLen)
			{
				ScoreType seedScore = oldSeedScore+1;
				result.VP |= (Word)1 << (Word)offset;
				oldSeedScore = seedScore;
				continue;
			}
			ScoreType seedScore = ((ScoreType)(j+offset+1)*100 + params.XscoreErrorCost / 2) / params.XscoreErrorCost;
			assert((size_t)seedScore < j + WordConfiguration<Word>::WordSize);
			assert((size_t)oldSeedScore < j + WordConfiguration<Word>::WordSize);
			assert(seedScore == oldSeedScore || seedScore == oldSeedScore+1);
			if (seedScore == oldSeedScore+1)
			{
				result.VP |= (Word)1 << (Word)offset;
			}
			oldSeedScore = seedScore;
		}
		result.scoreEnd = oldSeedScore;
		result.VP &= ~((Word)1);
		result.VN |= (Word)1;
		assert(result.getScoreBeforeStart() < (ScoreType)j + (ScoreType)WordConfiguration<Word>::WordSize);
		return result;
	}

	class StartComparer
	{
	public:
		StartComparer(const DPTable& slice) :
		slice(slice)
		{
		}
		bool operator()(const std::pair<size_t, size_t> left, const std::pair<size_t, size_t> right) const
		{
			return slice.slices[left.first].nodeMaxExactEndposScore.at(left.second) > slice.slices[right.first].nodeMaxExactEndposScore.at(right.second);
		}
	private:
		const DPTable& slice;
	};

	static std::vector<OnewayTrace> getLocalMaximaTracesFromTable(const Params& params, const std::string_view& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState, bool sliceConsistency, bool multiseed, std::vector<ScoreType>& sliceMaxScores)
	{
		std::vector<OnewayTrace> result;
		assert(slice.slices.size() > 1);
		assert(slice.slices[0].scores.size() == 0);
		std::vector<ScoreType> localSliceMaxScores;
		localSliceMaxScores.resize(slice.slices.size());
		for (size_t i = 0; i < slice.slices.size(); i++)
		{
			localSliceMaxScores[i] = slice.slices[i].maxExactEndposScore;
		}
		assert(localSliceMaxScores.size() <= sliceMaxScores.size());
		std::vector<std::pair<size_t, size_t>> possibleStarts;
		for (size_t currentSlice = slice.slices.size()-1; currentSlice > 0; currentSlice--)
		{
			assert(slice.slices[currentSlice].nodeMaxExactEndposScore.size() == slice.slices[currentSlice].scores.size());
			if (slice.slices[currentSlice].maxExactEndposScore <= 0) continue;
			std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, StartComparer> possibleStartsQueue { StartComparer { slice } };
			for (auto pair : slice.slices[currentSlice].scores)
			{
				if (!pair.second.exists) continue;
				auto node = pair.first;
				assert(slice.slices[currentSlice].nodeMaxExactEndposScore.count(node) == 1);
				ScoreType score = slice.slices[currentSlice].nodeMaxExactEndposScore.at(node);
				if (score <= 0) continue;
				if (score < localSliceMaxScores[currentSlice] * params.multimapScoreFraction) continue;
				if (score < sliceMaxScores[currentSlice] * params.multimapScoreFraction) continue;
				assert(currentSlice > 0);
				if (slice.slices[currentSlice-1].nodeMaxExactEndposScore.count(node) == 1)
				{
					if (slice.slices[currentSlice-1].nodeMaxExactEndposScore.at(node) > score) continue;
				}
				if (currentSlice < slice.slices.size()-1 && slice.slices[currentSlice+1].nodeMaxExactEndposScore.count(node) == 1)
				{
					if (slice.slices[currentSlice+1].nodeMaxExactEndposScore.at(node) >= score) continue;
				}
				bool notLocalMaximum = false;
				for (auto neighbor : params.graph.InNeighbors(node))
				{
					if (slice.slices[currentSlice].nodeMaxExactEndposScore.count(neighbor) == 1)
					{
						if (slice.slices[currentSlice].nodeMaxExactEndposScore.at(neighbor) > score)
						{
							notLocalMaximum = true;
							break;
						}
					}
					assert(currentSlice > 0);
					if (slice.slices[currentSlice-1].nodeMaxExactEndposScore.count(neighbor) == 1)
					{
						if (slice.slices[currentSlice-1].nodeMaxExactEndposScore.at(neighbor) > score)
						{
							notLocalMaximum = true;
							break;
						}
					}
				}
				for (auto neighbor : params.graph.OutNeighbors(node))
				{
					if (slice.slices[currentSlice].nodeMaxExactEndposScore.count(neighbor) == 1)
					{
						if (slice.slices[currentSlice].nodeMaxExactEndposScore.at(neighbor) > score)
						{
							notLocalMaximum = true;
							break;
						}
					}
					if (currentSlice < slice.slices.size()-1 && slice.slices[currentSlice+1].nodeMaxExactEndposScore.count(neighbor) == 1)
					{
						if (slice.slices[currentSlice+1].nodeMaxExactEndposScore.at(neighbor) > score)
						{
							notLocalMaximum = true;
							break;
						}
					}
				}
				if (notLocalMaximum) continue;
				if (possibleStartsQueue.size() < params.maxTraceCount)
				{
					possibleStartsQueue.emplace(currentSlice, node);
				}
				else
				{
					if (slice.slices[currentSlice].nodeMaxExactEndposScore.at(node) > slice.slices[possibleStartsQueue.top().first].nodeMaxExactEndposScore.at(possibleStartsQueue.top().second))
					{
						possibleStartsQueue.pop();
						possibleStartsQueue.emplace(currentSlice, node);
					}
				}
			}
			assert(possibleStartsQueue.size() <= params.maxTraceCount);
			while (possibleStartsQueue.size() > 0)
			{
				possibleStarts.emplace_back(possibleStartsQueue.top());
				possibleStartsQueue.pop();
			}
		}
		std::sort(possibleStarts.begin(), possibleStarts.end(), [&slice](const std::pair<size_t, size_t> left, const std::pair<size_t, size_t> right) { return slice.slices[left.first].nodeMaxExactEndposScore.at(left.second) > slice.slices[right.first].nodeMaxExactEndposScore.at(right.second); } );
		std::vector<OnewayTrace> fakeTraces;
		for (size_t starti = 0; starti < possibleStarts.size(); starti++)
		{
			assert(starti == 0 || slice.slices[possibleStarts[starti].first].nodeMaxExactEndposScore.at(possibleStarts[starti].second) <= slice.slices[possibleStarts[starti-1].first].nodeMaxExactEndposScore.at(possibleStarts[starti-1].second));
			size_t currentSlice = possibleStarts[starti].first;
			size_t node = possibleStarts[starti].second;
			ScoreType score = slice.slices[currentSlice].nodeMaxExactEndposScore.at(node);
			if (score < localSliceMaxScores[currentSlice] * params.multimapScoreFraction) continue;
			typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem previous;
			assert(currentSlice > 0);
			if (slice.slices[currentSlice-1].scores.hasNode(node))
			{
				previous = slice.slices[currentSlice-1].scores.node(node);
			}
			else
			{
				for (size_t i = 0; i < previous.NUM_CHUNKS; i++)
				{
					previous.HP[i] = WordConfiguration<Word>::AllOnes;
					previous.HN[i] = WordConfiguration<Word>::AllZeros;
				}
			}
			EqVector EqV = getEqVector(sequence, slice.slices[currentSlice].j);
			WordSlice fakeSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
			WordSlice seedstartSlice = getSeedSlice(slice.slices[currentSlice].j, sequence.size(), params);
			WordSlice extraSlice = slice.slices[currentSlice].seedstartNodes.count(node) == 1 ? seedstartSlice : fakeSlice;
			const typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem& nodeInfo = slice.slices[currentSlice].scores.node(node);
			std::vector<WordSlice> nodeSlices = recalcNodeWordslice(params, node, nodeInfo, EqV, previous, sliceConsistency, extraSlice, slice.slices[currentSlice].j);
			assert(nodeSlices[0].scoreEnd == nodeInfo.startSlice.scoreEnd);
			assert(nodeSlices[0].VP == nodeInfo.startSlice.VP);
			assert(nodeSlices[0].VN == nodeInfo.startSlice.VN);
			assert(nodeSlices.back().scoreEnd == nodeInfo.endSlice.scoreEnd);
			assert(nodeSlices.back().VP == nodeInfo.endSlice.VP);
			assert(nodeSlices.back().VN == nodeInfo.endSlice.VN);
			size_t nodeOffset = std::numeric_limits<size_t>::max();
			size_t bvOffset = std::numeric_limits<size_t>::max();
			for (size_t i = 0; i < nodeSlices.size(); i++)
			{
				auto maxScore = nodeSlices[i].maxXScoreFirstSlices(slice.slices[currentSlice].j, params.XscoreErrorCost, std::min((size_t)WordConfiguration<Word>::WordSize, (size_t)(sequence.size() - slice.slices[currentSlice].j)));
				assert(maxScore <= score);
				if (maxScore == score)
				{
					for (size_t off = WordConfiguration<Word>::WordSize-1; off < WordConfiguration<Word>::WordSize; off--)
					{
						if (slice.slices[currentSlice].j + off >= sequence.size()) continue;
						auto scoreHere = nodeSlices[i].getXScore(slice.slices[currentSlice].j, off, params.XscoreErrorCost);
						assert(scoreHere <= score);
						if (scoreHere == score)
						{
							if (nodeOffset == std::numeric_limits<size_t>::max() || off > bvOffset)
							{
								nodeOffset = i;
								bvOffset = off;
							}
						}
					}
				}
			}
			assert(nodeOffset != std::numeric_limits<size_t>::max());
			assert(bvOffset != std::numeric_limits<size_t>::max());
			assert(slice.slices[currentSlice].j + bvOffset < sequence.size());
			ScoreType startScore = nodeSlices[nodeOffset].getValue(bvOffset);
			MatrixPosition startPos { node, nodeOffset, slice.slices[currentSlice].j + bvOffset };
			auto oneTrace = getReverseTraceFromTable(params, result, fakeTraces, sequence, slice, reusableState, startPos, startScore, sliceConsistency, multiseed);
			if (!oneTrace.failed())
			{
				bool skip = false;
				for (size_t i = 0; i < result.size(); i++)
				{
					if (result[i].trace.back().DPposition != oneTrace.trace.back().DPposition) continue;
					if (result[i].alignmentXScore(params.XscoreErrorCost) >= oneTrace.alignmentXScore(params.XscoreErrorCost))
					{
						skip = true;
					}
					else
					{
						std::swap(result[i], result.back());
						result.pop_back();
						break;
					}
				}
				if (!skip)
				{
					result.emplace_back(std::move(oneTrace));
					size_t traceEnd = result.back().trace[0].DPposition.seqPos;
					size_t traceStart = result.back().trace.back().DPposition.seqPos;
					assert(traceEnd >= traceStart);
					size_t startSlice = ((traceEnd - traceStart)*0.05 + traceStart + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize + 1;
					for (size_t i = startSlice; i < currentSlice; i++)
					{
						localSliceMaxScores[i] = std::max(localSliceMaxScores[i], score);
					}
				}
			}
		}
		for (size_t i = 0; i < localSliceMaxScores.size(); i++)
		{
			sliceMaxScores[i] = std::max(sliceMaxScores[i], (ScoreType)(localSliceMaxScores[i] * 0.5));
		}
		return result;
	}

	static OnewayTrace getReverseTraceFromTableExactEndPos(const Params& params, const std::string_view& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState, bool sliceConsistency, bool multiseed)
	{
		assert(slice.slices.size() > 1);
		size_t bestIndex = 1;
		assert(slice.slices[1].maxExactEndposScore != std::numeric_limits<ScoreType>::max());
		for (size_t i = 1; i < slice.slices.size(); i++)
		{
			assert(slice.slices[i].maxExactEndposScore != std::numeric_limits<ScoreType>::max());
			if (slice.slices[i].maxExactEndposScore > slice.slices[bestIndex].maxExactEndposScore)
			{
				bestIndex = i;
			}
		}
		auto node = slice.slices[bestIndex].maxExactEndposNode;
		auto score = slice.slices[bestIndex].maxExactEndposScore;
		typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem previous;
		if (slice.slices[bestIndex-1].scores.hasNode(node))
		{
			previous = slice.slices[bestIndex-1].scores.node(node);
		}
		else
		{
			for (size_t i = 0; i < previous.NUM_CHUNKS; i++)
			{
				previous.HP[i] = WordConfiguration<Word>::AllOnes;
				previous.HN[i] = WordConfiguration<Word>::AllZeros;
			}
		}

		EqVector EqV = getEqVector(sequence, slice.slices[bestIndex].j);
		WordSlice fakeSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
		WordSlice seedstartSlice { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, std::numeric_limits<ScoreType>::max() };
		if (multiseed) seedstartSlice = getSeedSlice(slice.slices[bestIndex].j, sequence.size(), params);
		WordSlice extraSlice = slice.slices[bestIndex].seedstartNodes.count(node) == 1 ? seedstartSlice : fakeSlice;
		std::vector<WordSlice> nodeSlices = recalcNodeWordslice(params, node, slice.slices[bestIndex].scores.node(node), EqV, previous, sliceConsistency, extraSlice, slice.slices[bestIndex].j);

		size_t nodeOffset = std::numeric_limits<size_t>::max();
		size_t bvOffset = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < nodeSlices.size(); i++)
		{
			auto maxScore = nodeSlices[i].maxXScoreFirstSlices((ScoreType)slice.slices[bestIndex].j, params.XscoreErrorCost, std::min((size_t)WordConfiguration<Word>::WordSize, (size_t)(sequence.size() - slice.slices[bestIndex].j)));
			assert(maxScore <= score);
			if (maxScore == score)
			{
				for (size_t off = WordConfiguration<Word>::WordSize-1; off < WordConfiguration<Word>::WordSize; off--)
				{
					if (slice.slices[bestIndex].j + off >= sequence.size()) continue;
					auto scoreHere = nodeSlices[i].getXScore((ScoreType)slice.slices[bestIndex].j, off, params.XscoreErrorCost);
					assert(scoreHere <= score);
					if (scoreHere == score)
					{
						if (nodeOffset == std::numeric_limits<size_t>::max() || off > bvOffset)
						{
							nodeOffset = i;
							bvOffset = off;
						}
					}
				}
			}
		}
		assert(nodeOffset != std::numeric_limits<size_t>::max());
		assert(bvOffset != std::numeric_limits<size_t>::max());
		assert(slice.slices[bestIndex].j + bvOffset < sequence.size());
		ScoreType startScore = nodeSlices[nodeOffset].getValue(bvOffset);
		MatrixPosition startPos { node, nodeOffset, slice.slices[bestIndex].j + bvOffset };
		return getReverseTraceFromTable(params, sequence, slice, reusableState, startPos, startScore, sliceConsistency, multiseed);
	}

	static OnewayTrace getReverseTraceFromTableStartLastRow(const Params& params, const std::string_view& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState, bool sliceConsistency, bool multiseed)
	{
		ScoreType startScore = slice.slices.back().minScore;
		MatrixPosition startPos {slice.slices.back().minScoreNode, slice.slices.back().minScoreNodeOffset, std::min(slice.slices.back().j + WordConfiguration<Word>::WordSize - 1, sequence.size()-1)};
		return getReverseTraceFromTable(params, sequence, slice, reusableState, startPos, startScore, sliceConsistency, multiseed);
	}

	static OnewayTrace getReverseTraceFromTable(const Params& params, const std::string_view& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState, MatrixPosition startPos, ScoreType startScore, bool sliceConsistency, bool multiseed)
	{
		std::vector<OnewayTrace> empty;
		std::vector<OnewayTrace> fakes;
		return getReverseTraceFromTable(params, empty, fakes, sequence, slice, reusableState, startPos, startScore, sliceConsistency, multiseed);
	}

	static OnewayTrace getReverseTraceFromTable(const Params& params, const std::vector<OnewayTrace>& otherTraces, std::vector<OnewayTrace>& fakeTraces, const std::string_view& sequence, const DPTable& slice, AlignerGraphsizedState& reusableState, MatrixPosition startPos, ScoreType startScore, bool sliceConsistency, bool multiseed)
	{
		assert(slice.slices.size() > 0);
		size_t currentSlice = startPos.seqPos / WordConfiguration<Word>::WordSize + 1;
		assert(currentSlice < slice.slices.size());
		assert(slice.slices[currentSlice].minScoreNode != std::numeric_limits<LengthType>::max());
		assert(slice.slices[currentSlice].minScoreNodeOffset != std::numeric_limits<LengthType>::max());
		OnewayTrace result;
		result.score = startScore;
		result.trace.emplace_back(startPos, false, sequence, params.graph);
		LengthType currentNode = std::numeric_limits<LengthType>::max();
		std::vector<WordSlice> nodeSlices;
		EqVector EqV = getEqVector(sequence, slice.slices[currentSlice].j);
		WordSlice extraSlice;
		ScoreType lastScore = std::numeric_limits<ScoreType>::max()-1;
		std::vector<size_t> otherTracePoses;
		std::vector<size_t> fakeTracePoses;
		for (size_t i = 0; i < otherTraces.size(); i++)
		{
			otherTracePoses.push_back(0);
		}
		for (size_t i = 0; i < fakeTraces.size(); i++)
		{
			fakeTracePoses.push_back(0);
		}
		while (result.trace.back().DPposition.seqPos != (size_t)-1)
		{
			for (size_t i = 0; i < otherTraces.size(); i++)
			{
				if (otherTracePoses[i] >= otherTraces[i].trace.size()) continue;
				while (otherTracePoses[i] < otherTraces[i].trace.size() && otherTraces[i].trace[otherTracePoses[i]].DPposition.seqPos > result.trace.back().DPposition.seqPos) otherTracePoses[i] += 1;
				if (otherTracePoses[i] >= otherTraces[i].trace.size()) continue;
				if (otherTraces[i].trace[otherTracePoses[i]].DPposition == result.trace.back().DPposition)
				{
					fakeTraces.emplace_back(std::move(result));
					return OnewayTrace::TraceFailed();
				}
			}
			for (size_t i = 0; i < fakeTraces.size(); i++)
			{
				if (fakeTracePoses[i] >= fakeTraces[i].trace.size()) continue;
				while (fakeTracePoses[i] < fakeTraces[i].trace.size() && fakeTraces[i].trace[fakeTracePoses[i]].DPposition.seqPos > result.trace.back().DPposition.seqPos) fakeTracePoses[i] += 1;
				if (fakeTracePoses[i] >= fakeTraces[i].trace.size()) continue;
				if (fakeTraces[i].trace[fakeTracePoses[i]].DPposition == result.trace.back().DPposition)
				{
					fakeTraces.emplace_back(std::move(result));
					return OnewayTrace::TraceFailed();
				}
			}
			size_t newSlice = result.trace.back().DPposition.seqPos / WordConfiguration<Word>::WordSize + 1;
			assert(newSlice < slice.slices.size());
			assert(result.trace.back().DPposition.seqPos >= slice.slices[newSlice].j);
			assert(result.trace.back().DPposition.seqPos < slice.slices[newSlice].j + WordConfiguration<Word>::WordSize);
			LengthType newNode = result.trace.back().DPposition.node;
			if (newSlice != currentSlice || newNode != currentNode)
			{
				if (newSlice != currentSlice) EqV = getEqVector(sequence, slice.slices[newSlice].j);
				currentSlice = newSlice;
				currentNode = newNode;
				assert(slice.slices[currentSlice].scores.hasNode(currentNode));
				assert(currentSlice > 0);
				typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem previous;
				if (slice.slices[currentSlice-1].scores.hasNode(currentNode))
				{
					previous = slice.slices[currentSlice-1].scores.node(currentNode);
				}
				else
				{
					for (size_t i = 0; i < previous.NUM_CHUNKS; i++)
					{
						previous.HP[i] = WordConfiguration<Word>::AllOnes;
						previous.HN[i] = WordConfiguration<Word>::AllZeros;
					}
				}
				extraSlice.VP = WordConfiguration<Word>::AllZeros;
				extraSlice.VN = WordConfiguration<Word>::AllZeros;
				extraSlice.scoreEnd = std::numeric_limits<ScoreType>::max();
				if (slice.slices[currentSlice].seedstartNodes.count(currentNode) == 1) extraSlice = getSeedSlice(slice.slices[currentSlice].j, sequence.size(), params);
				nodeSlices = recalcNodeWordslice(params, currentNode, slice.slices[currentSlice].scores.node(currentNode), EqV, previous, sliceConsistency && !slice.slices[currentSlice].scoresNotValid && !slice.slices[currentSlice-1].scoresNotValid, extraSlice, slice.slices[currentSlice].j);
#ifdef SLICEVERBOSE
				std::cerr << "j " << slice.slices[currentSlice].j << " firstbt-calc " << slice.slices[currentSlice].scores.node(currentNode).firstSlicesCalcedWhenCalced << " lastbt-calc " << slice.slices[currentSlice].scores.node(currentNode).slicesCalcedWhenCalced << std::endl;
#endif
			}
			assert(result.trace.back().DPposition.node == currentNode);
			assert(result.trace.back().DPposition.nodeOffset < params.graph.NodeLength(currentNode));
			assert(nodeSlices.size() == params.graph.NodeLength(currentNode));
			assert(result.trace.back().DPposition.seqPos >= slice.slices[currentSlice].j);
			assert(result.trace.back().DPposition.seqPos < slice.slices[currentSlice].j + WordConfiguration<Word>::WordSize);
			assert((ScoreType)slice.slices[currentSlice].bandwidth >= 0);
			assert((ScoreType)slice.slices[currentSlice].bandwidth < std::numeric_limits<ScoreType>::max());
			assert(slice.slices[currentSlice].minScore < std::numeric_limits<ScoreType>::max() - (ScoreType)slice.slices[currentSlice].bandwidth);
			lastScore = nodeSlices[result.trace.back().DPposition.nodeOffset].getValue(result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize);
			if (result.trace.back().DPposition.seqPos == 0 && nodeSlices[result.trace.back().DPposition.nodeOffset].getValue(0) == extraSlice.getValue(0))
			{
				break;
			}
			if (result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize == 0 && result.trace.back().DPposition.nodeOffset == 0)
			{
				auto bt = pickBacktraceCorner(params, slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, currentNode, slice.slices[currentSlice].j, sequence, slice.slices[currentSlice].minScore + slice.slices[currentSlice].bandwidth, slice.slices[currentSlice].scoresNotValid, slice.slices[currentSlice-1].minScore + slice.slices[currentSlice-1].bandwidth, slice.slices[currentSlice-1].scoresNotValid, extraSlice);
				if (result.trace.size() > 0 && bt.first == result.trace.back().DPposition && !bt.second && nodeSlices[bt.first.nodeOffset].getValue(bt.first.seqPos % WordConfiguration<Word>::WordSize) == extraSlice.getValue(bt.first.seqPos % WordConfiguration<Word>::WordSize))
				{
					break;
				}
				if (result.trace.size() > 0 && bt.first == result.trace.back().DPposition && !bt.second && nodeSlices[bt.first.nodeOffset].getValue(bt.first.seqPos % WordConfiguration<Word>::WordSize) == extraSlice.getScoreBeforeStart())
				{
					break;
				}
				result.trace.emplace_back(bt.first, bt.second, sequence, params.graph);
				checkBacktraceCircularity(result);
				continue;
			}
			if (result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize == 0)
			{
				assert(currentSlice > 0);
				assert(result.trace.back().DPposition.nodeOffset > 0);
				if (!slice.slices[currentSlice-1].scores.hasNode(currentNode))
				{
					for (size_t i = result.trace.back().DPposition.nodeOffset-1; i > 0; i--)
					{
						result.trace.emplace_back(MatrixPosition {currentNode, i, result.trace.back().DPposition.seqPos}, false, sequence, params.graph);
					}
					result.trace.emplace_back(MatrixPosition {currentNode, 0, result.trace.back().DPposition.seqPos}, false, sequence, params.graph);
					continue;
				}
				auto crossing = pickBacktraceVerticalCrossing(params, slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, nodeSlices, slice.slices[currentSlice].j, currentNode, result.trace.back().DPposition, sequence, slice.slices[currentSlice].minScore + slice.slices[currentSlice].bandwidth, slice.slices[currentSlice].scoresNotValid, slice.slices[currentSlice-1].minScore + slice.slices[currentSlice-1].bandwidth, slice.slices[currentSlice-1].scoresNotValid, extraSlice);
				if (result.trace.size() > 0 && crossing.second.first.seqPos == result.trace.back().DPposition.seqPos && crossing.second.first.node == result.trace.back().DPposition.node && nodeSlices[crossing.second.first.nodeOffset].getValue(crossing.second.first.seqPos % WordConfiguration<Word>::WordSize) == extraSlice.getValue(crossing.second.first.seqPos % WordConfiguration<Word>::WordSize))
				{
					break;
				}
				assert(crossing.first.first.node == result.trace.back().DPposition.node);
				assert(crossing.first.first.seqPos == result.trace.back().DPposition.seqPos);
				assert(crossing.first.first.nodeOffset <= result.trace.back().DPposition.nodeOffset);
				assert(!crossing.first.second);
				if (crossing.first.first.nodeOffset != result.trace.back().DPposition.nodeOffset)
				{
					for (size_t nodeOffset = result.trace.back().DPposition.nodeOffset-1; nodeOffset != crossing.first.first.nodeOffset; nodeOffset--)
					{
						result.trace.emplace_back(MatrixPosition { crossing.first.first.node, nodeOffset, crossing.first.first.seqPos }, false, sequence, params.graph);
					}
				}
				if (crossing.first.first != result.trace.back().DPposition) result.trace.emplace_back(crossing.first.first, crossing.first.second, sequence, params.graph);
				assert(crossing.first.first == result.trace.back().DPposition);
				assert(crossing.second.first != result.trace.back().DPposition);
				result.trace.emplace_back(crossing.second.first, crossing.second.second, sequence, params.graph);
				continue;
			}
			if (result.trace.back().DPposition.nodeOffset == 0)
			{
				assert(result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize != 0);
				auto crossing = pickBacktraceHorizontalCrossing(params, slice.slices[currentSlice].scores, slice.slices[currentSlice-1].scores, slice.slices[currentSlice].j, currentNode, result.trace.back().DPposition, sequence, slice.slices[currentSlice].minScore + slice.slices[currentSlice].bandwidth, slice.slices[currentSlice].scoresNotValid, slice.slices[currentSlice-1].minScore + slice.slices[currentSlice-1].bandwidth, slice.slices[currentSlice-1].scoresNotValid, extraSlice);
				if (result.trace.size() > 0 && crossing.second.first.node == result.trace.back().DPposition.node && !crossing.second.second && nodeSlices[crossing.second.first.nodeOffset].getValue(crossing.second.first.seqPos % WordConfiguration<Word>::WordSize) == extraSlice.getValue(crossing.second.first.seqPos % WordConfiguration<Word>::WordSize))
				{
					break;
				}
				assert(crossing.first.first.node == result.trace.back().DPposition.node);
				assert(crossing.first.first.nodeOffset == result.trace.back().DPposition.nodeOffset);
				assert(crossing.first.first.seqPos <= result.trace.back().DPposition.seqPos);
				assert(!crossing.first.second);
				if (crossing.first.first.seqPos != result.trace.back().DPposition.seqPos)
				{
					for (size_t seqPos = result.trace.back().DPposition.seqPos-1; seqPos != crossing.first.first.seqPos; seqPos--)
					{
						result.trace.emplace_back(MatrixPosition { crossing.first.first.node, crossing.first.first.nodeOffset, seqPos }, false, sequence, params.graph);
					}
				}
				if (crossing.first.first != result.trace.back().DPposition) result.trace.emplace_back(crossing.first.first, crossing.first.second, sequence, params.graph);
				assert(crossing.first.first == result.trace.back().DPposition);
				assert(crossing.second.first != result.trace.back().DPposition);
				result.trace.emplace_back(crossing.second.first, crossing.second.second, sequence, params.graph);
				checkBacktraceCircularity(result);
				continue;
			}
			assert(result.trace.back().DPposition.nodeOffset != 0);
			assert(result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize != 0);
			auto inner = pickBacktraceInside(params, slice.slices[currentSlice].j, nodeSlices, result.trace.back().DPposition, sequence, extraSlice);
			for (auto pos : inner)
			{
				result.trace.emplace_back(pos, false, sequence, params.graph);
			}
			if (nodeSlices[result.trace.back().DPposition.nodeOffset].getValue(result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize) == extraSlice.getValue(result.trace.back().DPposition.seqPos % WordConfiguration<Word>::WordSize))
			{
				break;
			}
		}
		assert(result.score >= lastScore);
		result.score -= lastScore;
		do
		{
			if (multiseed) break;
			assert(result.trace.back().DPposition.seqPos == (size_t)-1);
			assert(slice.slices[0].scores.hasNode(result.trace.back().DPposition.node));
			auto node = slice.slices[0].scores.node(result.trace.back().DPposition.node);
			std::vector<ScoreType> beforeSliceScores;
			beforeSliceScores.resize(params.graph.NodeLength(result.trace.back().DPposition.node));
			beforeSliceScores[0] = node.startSlice.scoreEnd;
			for (size_t i = 1; i < beforeSliceScores.size(); i++)
			{
				size_t chunk = i / params.graph.SPLIT_NODE_SIZE;
				size_t offset = i % params.graph.SPLIT_NODE_SIZE;
				Word mask = ((Word)1) << offset;
				beforeSliceScores[i] = beforeSliceScores[i-1] + ((node.HP[chunk] & mask) >> offset) - ((node.HN[chunk] & mask) >> offset);
			}
			assert(beforeSliceScores.back() == node.endSlice.scoreEnd);
			while (beforeSliceScores[result.trace.back().DPposition.nodeOffset] != 0 && result.trace.back().DPposition.nodeOffset > 0 && beforeSliceScores[result.trace.back().DPposition.nodeOffset-1] == beforeSliceScores[result.trace.back().DPposition.nodeOffset] - 1)
			{
				result.trace.emplace_back(MatrixPosition {result.trace.back().DPposition.node, result.trace.back().DPposition.nodeOffset-1, result.trace.back().DPposition.seqPos}, false, sequence, params.graph);
			}
			if (result.trace.back().DPposition.nodeOffset == 0 && beforeSliceScores[result.trace.back().DPposition.nodeOffset] != 0)
			{
				bool found = false;
				for (auto neighbor : params.graph.InNeighbors(result.trace.back().DPposition.node))
				{
					if (slice.slices[0].scores.hasNode(neighbor) && slice.slices[0].scores.node(neighbor).endSlice.getScoreBeforeStart() == beforeSliceScores[result.trace.back().DPposition.nodeOffset] - 1)
					{
						result.trace.emplace_back(MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, result.trace.back().DPposition.seqPos}, true, sequence, params.graph);
						found = true;
						break;
					}
				}
				if (found) continue;
			}
		} while (false);
		return result;
	}

	static void checkBacktraceCircularity(const OnewayTrace& result)
	{
		for (size_t i = result.trace.size()-2; i < result.trace.size(); i--)
		{
			assert(result.trace[i].DPposition != result.trace.back().DPposition);
			if (result.trace[i].DPposition == result.trace.back().DPposition) std::abort();
			if (result.trace[i].DPposition.seqPos != result.trace.back().DPposition.seqPos) return;
		}
	}

	static std::vector<MatrixPosition> pickBacktraceInside(const Params& params, LengthType verticalOffset, const std::vector<WordSlice>& nodeSlices, MatrixPosition pos, const std::string_view& sequence, const WordSlice extraSlice)
	{
		assert(verticalOffset <= pos.seqPos);
		assert(verticalOffset + WordConfiguration<Word>::WordSize > pos.seqPos);
		assert((verticalOffset % WordConfiguration<Word>::WordSize) == 0);
		size_t hori = pos.nodeOffset;
		size_t vert = pos.seqPos - verticalOffset;
		assert(vert >= 0);
		assert(vert < WordConfiguration<Word>::WordSize);
		assert(hori >= 0);
		assert(hori < nodeSlices.size());
		std::vector<MatrixPosition> result;
		while (hori > 0 && vert > 0)
		{
			ScoreType scoreHere = nodeSlices[hori].getValue(vert);
			ScoreType verticalScore = nodeSlices[hori].getValue(vert-1);
			ScoreType horizontalScore = nodeSlices[hori-1].getValue(vert);
			ScoreType diagonalScore = nodeSlices[hori-1].getValue(vert-1);
			bool eq = Common::characterMatch(sequence[vert + verticalOffset], params.graph.NodeSequences(pos.node, hori));
			assert(verticalScore >= scoreHere-1);
			assert(horizontalScore >= scoreHere-1);
			assert(diagonalScore >= scoreHere - (eq?0:1));
			if (diagonalScore == scoreHere - (eq?0:1))
			{
				hori--;
				vert--;
				result.emplace_back(pos.node, hori, vert+verticalOffset);
				continue;
			}
			else if (verticalScore == scoreHere - 1)
			{
				vert--;
				result.emplace_back(pos.node, hori, vert+verticalOffset);
				continue;
			}
			else if (horizontalScore == scoreHere - 1)
			{
				hori--;
				result.emplace_back(pos.node, hori, vert+verticalOffset);
				continue;
			}
			assert(scoreHere == extraSlice.getValue(vert));
			break;
		}
		return result;
	}

	static std::pair<std::pair<MatrixPosition, bool>, std::pair<MatrixPosition, bool>> pickBacktraceHorizontalCrossing(const Params& params, const NodeSlice<LengthType, ScoreType, Word, false>& current, const NodeSlice<LengthType, ScoreType, Word, false>& previous, size_t j, LengthType node, MatrixPosition pos, const std::string_view& sequence, ScoreType quitScore, bool scoresNotValid, ScoreType previousQuitScore, bool previousScoresNotValid, const WordSlice extraSlice)
	{
		assert(current.hasNode(node));
		auto startSlice = current.node(node).startSlice;
		while (pos.seqPos % WordConfiguration<Word>::WordSize != 0 && (startSlice.VP & ((Word)1 << (pos.seqPos % WordConfiguration<Word>::WordSize))))
		{
			pos.seqPos--;
		}
		size_t offset = pos.seqPos % WordConfiguration<Word>::WordSize;
		if (offset == 0)
		{
			return std::make_pair(std::make_pair(pos, false), pickBacktraceCorner(params, current, previous, node, j, sequence, quitScore, scoresNotValid, previousQuitScore, previousScoresNotValid, extraSlice));
		}
		bool eq = Common::characterMatch(sequence[pos.seqPos], params.graph.NodeSequences(pos.node, pos.nodeOffset));
		ScoreType scoreHere = startSlice.getValue(offset);
		if (scoresNotValid || scoreHere > quitScore)
		{
			//this location is out of the band so the usual horizontal and vertical score limits don't apply
			//just pick the smallest scoring in-neighbor
			assert(offset > 0);
			ScoreType smallestFound = startSlice.getValue(offset-1);
			MatrixPosition smallestPos { node, 0, pos.seqPos-1 };
			bool nodeChange = false;
			for (auto neighbor : params.graph.InNeighbors(node))
			{
				if (current.hasNode(neighbor))
				{
					auto neighborSlice = current.node(neighbor).endSlice;
					if (neighborSlice.getValue(offset-1) <= smallestFound)
					{
						smallestFound = neighborSlice.getValue(offset-1);
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos-1 };
						nodeChange = true;
					}
					if (neighborSlice.getValue(offset) < smallestFound && neighbor != node)
					{
						smallestFound = neighborSlice.getValue(offset);
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos };
						nodeChange = true;
					}
				}
			}
			if (smallestPos == pos)
			{
				assert(extraSlice.getValue(offset) == scoreHere);
			}
			return std::make_pair(std::make_pair(pos, false), std::make_pair(smallestPos, nodeChange));
		}
		for (auto neighbor : params.graph.InNeighbors(node))
		{
			if (current.hasNode(neighbor))
			{
				auto neighborSlice = current.node(neighbor).endSlice;
				assert(neighborSlice.getValue(offset) >= scoreHere-1);
				assert(neighborSlice.getValue(offset-1) >= scoreHere - (eq?0:1));
				if (neighborSlice.getValue(offset-1) == scoreHere - (eq?0:1))
				{
					return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos-1 }, true));
				}
				if (neighborSlice.getValue(offset) == scoreHere-1)
				{
					return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, pos.seqPos }, true));
				}
			}
		}
		assert(scoreHere == extraSlice.getValue(offset));
		return std::make_pair(std::make_pair(pos, false), std::make_pair(pos, false));
	}

	static std::pair<std::pair<MatrixPosition, bool>, std::pair<MatrixPosition, bool>> pickBacktraceVerticalCrossing(const Params& params, const NodeSlice<LengthType, ScoreType, Word, false>& current, const NodeSlice<LengthType, ScoreType, Word, false>& previous, const std::vector<WordSlice> nodeScores, size_t j, LengthType node, MatrixPosition pos, const std::string_view& sequence, ScoreType quitScore, bool scoresNotValid, ScoreType previousQuitScore, bool previousScoresNotValid, const WordSlice extraSlice)
	{
		assert(pos.nodeOffset > 0);
		assert(pos.nodeOffset < nodeScores.size());
		while (pos.nodeOffset > 0 && nodeScores[pos.nodeOffset-1].getValue(0) == nodeScores[pos.nodeOffset].getValue(0) - 1)
		{
			pos.nodeOffset--;
		}
		if (pos.nodeOffset == 0)
		{
			return std::make_pair(std::make_pair(pos, false), pickBacktraceCorner(params, current, previous, node, j, sequence, quitScore, scoresNotValid, previousQuitScore, previousScoresNotValid, extraSlice));
		}
		assert(previous.hasNode(node));
		bool eq = Common::characterMatch(sequence[pos.seqPos], params.graph.NodeSequences(pos.node, pos.nodeOffset));
		auto previousNode = previous.node(node);
		ScoreType scoreHere = nodeScores[pos.nodeOffset].getValue(0);
		ScoreType scoreDiagonal = previousNode.startSlice.scoreEnd;
		for (size_t i = 1; i <= pos.nodeOffset - 1; i++)
		{
			scoreDiagonal += (previousNode.HP[i / WordConfiguration<Word>::WordSize] >> (i % WordConfiguration<Word>::WordSize)) & 1;
			scoreDiagonal -= (previousNode.HN[i / WordConfiguration<Word>::WordSize] >> (i % WordConfiguration<Word>::WordSize)) & 1;
		}
		ScoreType scoreUp = scoreDiagonal;
		scoreUp += (previousNode.HP[(pos.nodeOffset) / WordConfiguration<Word>::WordSize] >> ((pos.nodeOffset) % WordConfiguration<Word>::WordSize)) & 1;
		scoreUp -= (previousNode.HN[(pos.nodeOffset) / WordConfiguration<Word>::WordSize] >> ((pos.nodeOffset) % WordConfiguration<Word>::WordSize)) & 1;
		if (previousScoresNotValid || scoresNotValid || scoreHere > quitScore || scoreDiagonal > previousQuitScore || scoreUp > previousQuitScore)
		{
			if (scoreHere == extraSlice.getValue(0))
			{
				return std::make_pair(std::make_pair(pos, false), std::make_pair(pos, false));
			}
			//this location is out of the band so the usual horizontal and vertical score limits don't apply
			//just pick the smallest scoring in-neighbor
			if (scoreDiagonal <= scoreUp)
			{
				return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset - 1, pos.seqPos-1}, false));
			}
			else
			{
				return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset, pos.seqPos-1}, false));
			}
		}
		assert(scoreUp >= scoreHere - 1);
		assert(scoreDiagonal >= scoreHere - (eq?0:1));
		if (scoreDiagonal == scoreHere - (eq?0:1)) return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset - 1, pos.seqPos-1}, false));
		if (scoreUp == scoreHere - 1) return std::make_pair(std::make_pair(pos, false), std::make_pair(MatrixPosition{pos.node, pos.nodeOffset, pos.seqPos-1}, false));
		assert(scoreHere == extraSlice.getValue(0));
		return std::make_pair(std::make_pair(pos, false), std::make_pair(pos, false));
	}

	static std::pair<MatrixPosition, bool> pickBacktraceCorner(const Params& params, const NodeSlice<LengthType, ScoreType, Word, false>& current, const NodeSlice<LengthType, ScoreType, Word, false>& previous, LengthType node, size_t j, const std::string_view& sequence, ScoreType quitScore, bool scoresNotValid, ScoreType previousQuitScore, bool previousScoresNotValid, const WordSlice extraSlice)
	{
		ScoreType scoreHere = current.node(node).startSlice.getValue(0);
		if (scoresNotValid || scoreHere > quitScore)
		{
			if (scoreHere == extraSlice.getValue(0))
			{
				return std::make_pair(MatrixPosition { node, 0, j }, false);
			}
			//this location is out of the band so the usual horizontal and vertical score limits don't apply
			//just pick the smallest scoring in-neighbor
			ScoreType smallestFound = scoreHere+1;
			MatrixPosition smallestPos { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
			bool nodeChange = false;
			if (previous.hasNode(node))
			{
				auto slice = previous.node(node).startSlice;
				smallestFound = slice.scoreEnd;
				smallestPos = MatrixPosition { node, 0, j-1 };
				nodeChange = false;
			}
			for (auto neighbor : params.graph.InNeighbors(node))
			{
				if (previous.hasNode(neighbor))
				{
					auto neighborSlice = previous.node(neighbor).endSlice;
					if (neighborSlice.scoreEnd <= smallestFound)
					{
						smallestFound = neighborSlice.scoreEnd;
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, j-1 };
						nodeChange = true;
					}
				}
				if (current.hasNode(neighbor) && neighbor != node)
				{
					auto neighborSlice = current.node(neighbor).endSlice;
					if (neighborSlice.getValue(0) < smallestFound)
					{
						smallestFound = neighborSlice.getValue(0);
						smallestPos = MatrixPosition { neighbor, params.graph.NodeLength(neighbor)-1, j };
						nodeChange = true;
					}
				}
			}
			assert(smallestPos.node != std::numeric_limits<size_t>::max());
			return std::make_pair(smallestPos, nodeChange);
		}
		assert(scoreHere <= quitScore);
		bool eq = Common::characterMatch(sequence[j], params.graph.NodeSequences(node, 0));
		MatrixPosition bestInvalidBacktrace { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
		ScoreType bestInvalidBacktraceScore = scoreHere+1;
		bool bestInvalidNodeSwitch = false;
		for (auto neighbor : params.graph.InNeighbors(node))
		{
			if (previous.hasNode(neighbor))
			{
				ScoreType cornerScore = previous.node(neighbor).endSlice.scoreEnd;
				if (previousScoresNotValid || cornerScore > previousQuitScore)
				{
					//scores not valid, pick best
					if (cornerScore < bestInvalidBacktraceScore)
					{
						bestInvalidBacktraceScore = cornerScore;
						bestInvalidBacktrace = MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, j-1 };
						bestInvalidNodeSwitch = true;
					}
				}
				else
				{
					assert(previous.node(neighbor).endSlice.scoreEnd >= scoreHere-(eq?0:1));
					if (cornerScore == scoreHere-(eq?0:1))
					{
						return std::make_pair(MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, j-1 }, true);
					}
				}
			}
			if (current.hasNode(neighbor))
			{
				assert(current.node(neighbor).endSlice.getValue(0) >= scoreHere-1);
				if (current.node(neighbor).endSlice.getValue(0) == scoreHere-1)
				{
					return std::make_pair(MatrixPosition {neighbor, params.graph.NodeLength(neighbor)-1, j}, true);
				}
			}
		}
		if (previous.hasNode(node))
		{
			assert(previous.node(node).startSlice.scoreEnd >= scoreHere-1);
			if (previous.node(node).startSlice.scoreEnd == scoreHere-1)
			{
				return std::make_pair(MatrixPosition {node, 0, j-1 }, false);
			}
			if (previousScoresNotValid || previous.node(node).startSlice.scoreEnd > previousQuitScore)
			{
				if (previous.node(node).startSlice.scoreEnd < bestInvalidBacktraceScore)
				{
					bestInvalidBacktraceScore = previous.node(node).startSlice.scoreEnd;
					bestInvalidBacktrace = MatrixPosition { node, 0, j-1 };
					bestInvalidNodeSwitch = false;
				}
			}
		}
		if (eq && scoreHere == extraSlice.getScoreBeforeStart())
		{
			return std::make_pair(MatrixPosition { node, 0, j }, false);
		}
		if (scoreHere == extraSlice.getValue(0))
		{
			return std::make_pair(MatrixPosition { node, 0, j }, false);
		}
		//scores not valid, pick best
		if (bestInvalidBacktraceScore < scoreHere+1)
		{
			assert(bestInvalidBacktrace.node != (size_t)-1 || bestInvalidBacktrace.nodeOffset != (size_t)-1 || bestInvalidBacktrace.seqPos != (size_t)-1);
			return std::make_pair(bestInvalidBacktrace, bestInvalidNodeSwitch);
		}
		assert(false);
		return std::make_pair(MatrixPosition {0, 0, 0}, false);
	}

	static WordSlice getSourceSliceFromScore(ScoreType previousScore)
	{
		WordSlice result { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, previousScore+WordConfiguration<Word>::WordSize };
		return result;
	}

	static void assertSliceCorrectness(WordSlice oldSlice, WordSlice newSlice, Word Eq, int hin, WordSlice extraSlice)
	{
		ScoreType foundMinScore = newSlice.getScoreBeforeStart() + 1;
		foundMinScore = std::min(foundMinScore, oldSlice.getValue(0) + 1);
		foundMinScore = std::min(foundMinScore, oldSlice.getScoreBeforeStart() + ((Eq & 1) ? 0 : 1));
		assert(newSlice.getScoreBeforeStart() == oldSlice.getScoreBeforeStart() + hin || newSlice.getScoreBeforeStart() == extraSlice.getScoreBeforeStart());
		assert(newSlice.getValue(0) == foundMinScore || newSlice.getValue(0) == extraSlice.getValue(0));
		for (size_t i = 1; i < WordConfiguration<Word>::WordSize; i++)
		{
			foundMinScore = newSlice.getValue(i-1)+1;
			foundMinScore = std::min(foundMinScore, oldSlice.getValue(i)+1);
			foundMinScore = std::min(foundMinScore, oldSlice.getValue(i-1) + ((Eq & ((Word)1 << i)) ? 0 : 1));
			assert(newSlice.getValue(i) == foundMinScore || newSlice.getValue(i) == extraSlice.getValue(i));
		}
	}

	static std::vector<WordSlice> recalcNodeWordslice(const Params& params, LengthType node, const typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem& slice, const EqVector& EqV, const typename NodeSlice<LengthType, ScoreType, Word, false>::NodeSliceMapItem& previousSlice, bool sliceConsistency, const WordSlice extraSlice, ScoreType seqOffset)
	{
		size_t nodeLength = params.graph.NodeLength(node);
		std::vector<EdgeWithPriority> incoming; // also fake!
		incoming.emplace_back(node, 0, slice.startSlice, true);
		std::vector<WordSlice> result;
		result.reserve(nodeLength);
		auto sliceCopy = slice;
		if (node < params.graph.FirstAmbiguous())
		{
			calculateNodeInner<false>(params, node, sliceCopy, EqV, previousSlice, incoming, [](size_t pos) { return false; }, params.graph.NodeChunks(node), extraSlice, [&result](const WordSlice& slice) { result.push_back(slice); }, seqOffset);
		}
		else
		{
			calculateNodeInner<false>(params, node, sliceCopy, EqV, previousSlice, incoming, [](size_t pos) { return false; }, params.graph.AmbiguousNodeChunks(node), extraSlice, [&result](const WordSlice& slice) { result.push_back(slice); }, seqOffset);
		}
		assert(result.size() == nodeLength);
		assert(!sliceConsistency || result[0].scoreEnd == slice.startSlice.scoreEnd);
		assert(!sliceConsistency || result[0].VP == slice.startSlice.VP);
		assert(!sliceConsistency || result[0].VN == slice.startSlice.VN);
		assert(!sliceConsistency || result.back().scoreEnd == slice.endSlice.scoreEnd);
		assert(!sliceConsistency || result.back().VP == slice.endSlice.VP);
		assert(!sliceConsistency || result.back().VN == slice.endSlice.VN);
		return result;
	}

	template <typename NodeChunkType>
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static NodeCalculationResult calculateNodeClipPrecise(const Params& params, size_t i, typename NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem& slice, const EqVector& EqV, typename NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem previousSlice, const std::vector<EdgeWithPriority>& incoming, const std::vector<bool>& previousBand, NodeChunkType nodeChunks, const WordSlice extraSlice, ScoreType seqOffset)
	{
		return calculateNodeInner<true>(params, i, slice, EqV, previousSlice, incoming, [&previousBand](size_t pos) { return previousBand[pos]; }, nodeChunks, extraSlice, [](const WordSlice& slice){}, seqOffset);
	}

	template <bool AllowEarlyLeave, typename NodeChunkType, typename WordsliceCallback, typename ExistenceCheckFunction>
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static NodeCalculationResult calculateNodeInner(const Params& params, size_t i, typename NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem& slice, const EqVector& EqV, typename NodeSlice<LengthType, ScoreType, Word, true>::NodeSliceMapItem previousSlice, const std::vector<EdgeWithPriority>& incoming, ExistenceCheckFunction bandCheck, NodeChunkType nodeChunks, const WordSlice extraSlice, WordsliceCallback callback, ScoreType seqOffset)
	{
		assert(incoming.size() > 0);
		WordSlice newWs;
		WordSlice ws;
		bool hasWs = false;
		NodeCalculationResult result;
		result.minScore = std::numeric_limits<ScoreType>::max();
		result.minScoreNode = std::numeric_limits<LengthType>::max();
		result.minScoreNodeOffset = std::numeric_limits<LengthType>::max();
		result.maxExactEndposScore = std::numeric_limits<ScoreType>::min();
		result.maxExactEndposNode = std::numeric_limits<LengthType>::max();
		result.cellsProcessed = 0;
		auto nodeLength = params.graph.NodeLength(i);

		Word Eq = EqV.getEqI(nodeChunks[0] & 3);
		bool hasSkipless = false;
		bool forceCalculation = false;

		for (auto inc : incoming)
		{
			result.cellsProcessed++;
			if (inc.forceCalculation) forceCalculation = true;
			if (inc.skipFirst)
			{
				if (!hasWs)
				{
					ws = inc.incoming;
					hasWs = true;
				}
				else
				{
					ws = ws.mergeWith(inc.incoming);
				}
				continue;
			}
			hasSkipless = true;
			Word hinP;
			Word hinN;
			ScoreType compareHin = std::numeric_limits<ScoreType>::max();
			if (previousSlice.exists)
			{
				compareHin = previousSlice.startSlice.scoreEnd;
			}
			if (extraSlice.scoreEnd != std::numeric_limits<ScoreType>::max())
			{
				compareHin = std::min(compareHin, extraSlice.getScoreBeforeStart());
			}
			ScoreType incomingScoreBeforeStart = inc.incoming.getScoreBeforeStart();
			if (compareHin < incomingScoreBeforeStart)
			{
				hinP = 0;
				hinN = 1;
			}
			else if (compareHin > incomingScoreBeforeStart)
			{
				hinP = 1;
				hinN = 0;
			}
			else
			{
				hinP = 0;
				hinN = 0;
			}

			WordSlice newWs;
			std::tie(newWs, hinP, hinN) = getNextSlice(Eq, inc.incoming, hinP, hinN);
			if (!previousSlice.exists || newWs.getScoreBeforeStart() < previousSlice.startSlice.scoreEnd)
			{
				newWs.VP &= WordConfiguration<Word>::AllOnes ^ 1;
				newWs.VN |= 1;
			}
			// assert(newWs.getScoreBeforeStart() >= debugLastRowMinScore || newWs.getScoreBeforeStart() >= extraSlice.getScoreBeforeStart());
			if (!hasWs)
			{
				ws = newWs;
				hasWs = true;
			}
			else
			{
				ws = ws.mergeWith(newWs);
			}
		}

		assert(hasWs);

		if (extraSlice.scoreEnd != std::numeric_limits<ScoreType>::max())
		{
			ws = ws.mergeWith(extraSlice);
		}

		result.minScore = ws.scoreEnd;
		result.minScoreNode = i;
		result.minScoreNodeOffset = 0;
		
		result.maxExactEndposScore = ws.maxXScore(seqOffset, params.XscoreErrorCost);
		result.maxExactEndposNode = i;

		if (!forceCalculation && slice.exists)
		{
			if (hasSkipless && params.graph.InNeighbors(i).size() == 1 && bandCheck(params.graph.InNeighbors(i)[0]))
			{
				if (ws.scoreEnd > slice.startSlice.scoreEnd)
				{
#ifdef EXTRACORRECTNESSASSERTIONS
					auto debugTest = ws.mergeWith(slice.startSlice);
					assert(debugTest.VP == slice.startSlice.VP);
					assert(debugTest.VN == slice.startSlice.VN);
					assert(debugTest.scoreEnd == slice.startSlice.scoreEnd);
#endif
					if constexpr (AllowEarlyLeave) return result;
				}
				else if (ws.scoreEnd < slice.startSlice.scoreEnd)
				{
#ifdef EXTRACORRECTNESSASSERTIONS
					auto debugTest = ws.mergeWith(slice.startSlice);
					// todo figure out why breaks
					assert(debugTest.VP == ws.VP);
					assert(debugTest.VN == ws.VN);
					assert(debugTest.scoreEnd == ws.scoreEnd);
#endif
				}
				else
				{
					Word newBigger = (ws.VP & ~slice.startSlice.VP) | (slice.startSlice.VN & ~ws.VN);
					Word oldBigger = (slice.startSlice.VP & ~ws.VP) | (ws.VN & ~slice.startSlice.VN);
					if (newBigger > oldBigger)
					{
#ifdef EXTRACORRECTNESSASSERTIONS
						auto debugTest = ws.mergeWith(slice.startSlice);
						// todo figure out why breaks
						assert(debugTest.VP == ws.VP);
						assert(debugTest.VN == ws.VN);
						assert(debugTest.scoreEnd == ws.scoreEnd);
#endif
					}
					else if (oldBigger > newBigger)
					{
#ifdef EXTRACORRECTNESSASSERTIONS
						auto debugTest = ws.mergeWith(slice.startSlice);
						assert(debugTest.VP == slice.startSlice.VP);
						assert(debugTest.VN == slice.startSlice.VN);
						assert(debugTest.scoreEnd == slice.startSlice.scoreEnd);
#endif
						if constexpr (AllowEarlyLeave) return result;
					}
					else if (newBigger == 0 && oldBigger == 0)
					{
						assert(ws.VP == slice.startSlice.VP);
						assert(ws.VN == slice.startSlice.VN);
						assert(ws.scoreEnd == slice.startSlice.scoreEnd);
						if constexpr (AllowEarlyLeave) return result;
					}
					else
					{
						WordSlice test = ws.mergeWith(slice.startSlice);
						if (test.scoreEnd == slice.startSlice.scoreEnd && test.VP == slice.startSlice.VP && test.VN == slice.startSlice.VN)
						{
							if constexpr (AllowEarlyLeave) return result;
						}
						ws = test;
					}
				}
			}
			else
			{
				WordSlice test = ws.mergeWith(slice.startSlice);
				if (test.scoreEnd == slice.startSlice.scoreEnd && test.VP == slice.startSlice.VP && test.VP == slice.startSlice.VN)
				{
					if constexpr (AllowEarlyLeave) return result;
				}
				ws = test;
			}
		}

		if (previousSlice.exists)
		{
			if (ws.getScoreBeforeStart() > previousSlice.startSlice.scoreEnd)
			{
				ws = ws.mergeWith(getSourceSliceFromScore(previousSlice.startSlice.scoreEnd));
			}
		}

		for (size_t i = 0; i < slice.NUM_CHUNKS; i++)
		{
			slice.HP[i] = WordConfiguration<Word>::AllZeros;
			slice.HN[i] = WordConfiguration<Word>::AllZeros;
		}

		LengthType pos = 1;
		Word forceMask = WordConfiguration<Word>::AllZeros;
		Word fixedHP[1] { previousSlice.HP[0] };
		Word fixedHN[1] { previousSlice.HN[0] };

		if (previousSlice.exists)
		{
			ScoreType scoreBefore = ws.getScoreBeforeStart();
			ScoreType scoreComparison = previousSlice.startSlice.scoreEnd;
			assert(scoreBefore <= scoreComparison);
			if (extraSlice.scoreEnd != std::numeric_limits<ScoreType>::max())
			{
				scoreBefore = std::min(scoreBefore, extraSlice.getScoreBeforeStart());
				assert(scoreBefore <= extraSlice.getScoreBeforeStart());
				assert(scoreBefore <= ws.getScoreBeforeStart());
				size_t fixoffset = 1;
				if (scoreBefore < scoreComparison) forceMask |= (Word)1;
				ScoreType adjustScore = std::min(scoreBefore, scoreComparison);
				for (size_t fixchunk = 0; fixchunk < slice.NUM_CHUNKS; fixchunk++)
				{
					for (; fixoffset < WordConfiguration<Word>::WordSize && fixchunk * WordConfiguration<Word>::WordSize + fixoffset < nodeLength; fixoffset++)
					{
						ScoreType newScoreComparison = scoreComparison;
						newScoreComparison += (fixedHP[fixchunk] >> fixoffset) & 1;
						newScoreComparison -= (fixedHN[fixchunk] >> fixoffset) & 1;
						if (scoreBefore < extraSlice.getScoreBeforeStart()) scoreBefore += 1;
						assert(scoreBefore <= extraSlice.getScoreBeforeStart());
						Word mask = ((Word)1) << fixoffset;
						if (scoreBefore < newScoreComparison)
						{
							forceMask |= mask;
						}
						ScoreType newAdjustScore = std::min(newScoreComparison, scoreBefore);
						assert(newAdjustScore >= adjustScore-1);
						assert(newAdjustScore <= adjustScore+1);
						if (newAdjustScore == adjustScore-1)
						{
							fixedHN[fixchunk] |= mask;
							fixedHP[fixchunk] &= ~mask;
						}
						else if (newAdjustScore == adjustScore)
						{
							fixedHP[fixchunk] &= ~mask;
							fixedHN[fixchunk] &= ~mask;
						}
						else if (newAdjustScore == adjustScore+1)
						{
							fixedHP[fixchunk] |= mask;
							fixedHN[fixchunk] &= ~mask;
						}
						adjustScore = newAdjustScore;
						scoreComparison = newScoreComparison;
					}
				}
				assert(scoreComparison == previousSlice.endSlice.scoreEnd);
				assert(adjustScore <= extraSlice.getScoreBeforeStart());
				assert(adjustScore <= previousSlice.endSlice.scoreEnd);
				assert(adjustScore <= ws.getScoreBeforeStart() + (ScoreType)nodeLength - 1);
				assert(adjustScore == previousSlice.endSlice.scoreEnd || adjustScore == extraSlice.getScoreBeforeStart() || adjustScore == ws.getScoreBeforeStart() + (ScoreType)nodeLength - 1);
			}
			else if (scoreBefore < scoreComparison)
			{
				size_t fixoffset = 1;
				forceMask |= (Word)1;
				for (size_t fixchunk = 0; fixchunk < slice.NUM_CHUNKS; fixchunk++)
				{
					for (; fixoffset < WordConfiguration<Word>::WordSize; fixoffset++)
					{
						ScoreType newScoreComparison = scoreComparison;
						newScoreComparison += (fixedHP[fixchunk] >> fixoffset) & 1;
						newScoreComparison -= (fixedHN[fixchunk] >> fixoffset) & 1;
						Word mask = ((Word)1) << fixoffset;
						assert(scoreBefore <= newScoreComparison);
						if (scoreBefore < newScoreComparison)
						{
							fixedHP[fixchunk] |= mask;
							fixedHN[fixchunk] &= ~mask;
							forceMask |= mask;
						}
						if (scoreBefore == newScoreComparison)
						{
							fixedHP[fixchunk] &= ~mask;
							fixedHN[fixchunk] &= ~mask;
						}
						scoreBefore++;
						scoreComparison = newScoreComparison;
						if (scoreBefore >= scoreComparison) break;
					}
					if (scoreBefore >= scoreComparison) break;
					fixoffset = 0;
				}
			}
		}
		else
		{
			forceMask = WordConfiguration<Word>::AllOnes;
		}
		if (ws.scoreEnd < result.minScore)
		{
			result.minScore = ws.scoreEnd;
			result.minScoreNodeOffset = pos;
		}
		bool endIncomparable = forceMask & ((Word)1 << (Word)(nodeLength-1));
		slice.startSlice = ws;
		if constexpr (!AllowEarlyLeave) callback(ws);
		slice.exists = true;
		Word forceEq = WordConfiguration<Word>::AllOnes;
		Word hinP, hinN;
		if (!previousSlice.exists) forceEq ^= 1;
		size_t smallChunk = 0;
		size_t offset = 1;
		pos = smallChunk * (WordConfiguration<Word>::WordSize / 2) + offset;
		for (; smallChunk < params.graph.CHUNKS_IN_NODE; smallChunk++)
		{
			size_t bigChunk = smallChunk / 2;
			size_t bigChunkOffset = (smallChunk % 2) * (WordConfiguration<Word>::WordSize / 2);
			Word HP = fixedHP[bigChunk] >> bigChunkOffset;
			Word HN = fixedHN[bigChunk] >> bigChunkOffset;
			auto charChunk = nodeChunks[smallChunk];
			HP >>= offset;
			HN >>= offset;
			charChunk >>= offset * 2;
			forceMask >>= offset;
			for (; offset < WordConfiguration<Word>::WordSize / 2 && pos < nodeLength; offset++)
			{
				Eq = EqV.getEqI(charChunk & 3);
				Eq &= forceEq;
				Word newHP;
				Word newHN;
				std::tie(newWs, newHP, newHN) = getNextSlice(Eq, ws, HP & 1, HN & 1);
				if (forceMask & 1)
				{
					newWs.VP &= WordConfiguration<Word>::AllOnes ^ 1;
					newWs.VN |= 1;
				}
				if (extraSlice.scoreEnd != std::numeric_limits<ScoreType>::max())
				{
					newWs = newWs.mergeWith(extraSlice);
					assert(newWs.scoreEnd >= ws.scoreEnd-1);
					assert(newWs.scoreEnd <= ws.scoreEnd+1);
					hinP = newWs.scoreEnd == ws.scoreEnd+1 ? 1 : 0;
					hinN = newWs.scoreEnd == ws.scoreEnd-1 ? 1 : 0;
				}
				else
				{
					hinP = newHP >> (WordConfiguration<Word>::WordSize-1);
					hinN = newHN >> (WordConfiguration<Word>::WordSize-1);
				}
#ifdef EXTRACORRECTNESSASSERTIONS
				assertSliceCorrectness(ws, newWs, Eq, (HP & 1) - (HN & 1), extraSlice);
				if (hinP == 1) assert(newWs.scoreEnd == ws.scoreEnd+1);
				if (hinN == 1) assert(newWs.scoreEnd == ws.scoreEnd-1);
#endif
				// assert(!AllowEarlyLeave || newWs.getScoreBeforeStart() >= debugLastRowMinScore || newWs.getScoreBeforeStart() >= extraSlice.getScoreBeforeStart());
				ws = newWs;
				if (ws.scoreEnd < result.minScore)
				{
					result.minScore = ws.scoreEnd;
					result.minScoreNodeOffset = pos;
				}
				result.maxExactEndposScore = std::max(result.maxExactEndposScore, ws.maxXScore(seqOffset, params.XscoreErrorCost, newHN));
				if constexpr (!AllowEarlyLeave) callback(ws);
				charChunk >>= 2;
				HP >>= 1;
				HN >>= 1;
				forceMask >>= 1;
				pos++;
				slice.HP[bigChunk] |= hinP << (offset + bigChunkOffset);
				slice.HN[bigChunk] |= hinN << (offset + bigChunkOffset);
			}
			offset = 0;
		}
		result.cellsProcessed = pos;
		slice.endSlice = ws;
#ifndef NDEBUG
		if (previousSlice.exists) assert(slice.endSlice.getScoreBeforeStart() <= previousSlice.endSlice.scoreEnd);
		if (previousSlice.exists && !endIncomparable) assert(slice.endSlice.getScoreBeforeStart() == previousSlice.endSlice.scoreEnd || slice.endSlice.getScoreBeforeStart() == extraSlice.getScoreBeforeStart());
#endif
		return result;
	}

	template <bool HasVectorMap, bool PreviousHasVectorMap>
	static void flattenLastSliceEnd(const Params& params, NodeSlice<LengthType, ScoreType, Word, HasVectorMap>& slice, const NodeSlice<LengthType, ScoreType, Word, PreviousHasVectorMap>& previousSlice, NodeCalculationResult& sliceCalc, LengthType j, const std::string_view& sequence, bool sliceConsistency, const WordSlice extraSlice)
	{
		assert(j < sequence.size());
		assert(sequence.size() - j < WordConfiguration<Word>::WordSize);
		sliceCalc.minScore = std::numeric_limits<ScoreType>::max();
		sliceCalc.minScoreNode = std::numeric_limits<LengthType>::max();
		sliceCalc.minScoreNodeOffset = std::numeric_limits<LengthType>::max();
		sliceCalc.maxExactEndposScore = std::numeric_limits<ScoreType>::min();
		sliceCalc.maxExactEndposNode = std::numeric_limits<LengthType>::max();
		auto offset = sequence.size() - j;
		assert(offset >= 0);
		assert(offset < WordConfiguration<Word>::WordSize);
		EqVector EqV = getEqVector(sequence, j);
		for (auto node : slice)
		{
			auto current = node.second;
			decltype(current) old;
			if (previousSlice.hasNode(node.first))
			{
				old = previousSlice.node(node.first);
			}
			else
			{
				old.exists = false;
				for (size_t i = 0; i < old.NUM_CHUNKS; i++)
				{
					old.HP[i] = WordConfiguration<Word>::AllOnes;
					old.HN[i] = WordConfiguration<Word>::AllZeros;
				}
			}

			std::vector<WordSlice> nodeSlices = recalcNodeWordslice(params, node.first, current, EqV, old, sliceConsistency, extraSlice, j);

			assert(nodeSlices[0].VP == node.second.startSlice.VP);
			assert(nodeSlices[0].VN == node.second.startSlice.VN);
			assert(nodeSlices[0].scoreEnd == node.second.startSlice.scoreEnd);
			assert(nodeSlices.back().VP == node.second.endSlice.VP);
			assert(nodeSlices.back().VN == node.second.endSlice.VN);
			assert(nodeSlices.back().scoreEnd == node.second.endSlice.scoreEnd);
			for (size_t i = 0; i < nodeSlices.size(); i++)
			{
				auto wordSliceResult = flattenWordSlice(nodeSlices[i], offset);
				if (wordSliceResult.scoreEnd < sliceCalc.minScore)
				{
					sliceCalc.minScore = wordSliceResult.scoreEnd;
					sliceCalc.minScoreNode = node.first;
					sliceCalc.minScoreNodeOffset = i;
				}
				if (wordSliceResult.maxXScoreFirstSlices(j, params.XscoreErrorCost, offset) > sliceCalc.maxExactEndposScore)
				{
					sliceCalc.maxExactEndposScore = wordSliceResult.maxXScore(j, params.XscoreErrorCost);
					sliceCalc.maxExactEndposNode = i;
				}
			}
		}
		assert(sliceCalc.minScore != std::numeric_limits<ScoreType>::max());
		assert(sliceCalc.minScoreNode < params.graph.NodeSize());
		assert(sliceCalc.minScoreNodeOffset < params.graph.NodeLength(sliceCalc.minScoreNode));
	}

	static DPSlice getInitialEmptySlice()
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		result.scores.addEmptyNodeMap(1);
		return result;
	}

	static DPSlice getInitialSliceExactPosition(const Params& params, LengthType bigraphNodeId, size_t offset)
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		result.scores.addEmptyNodeMap(1);
		assert(offset < params.graph.BigraphNodeSize(bigraphNodeId));
		size_t nodeIndex = params.graph.GetDigraphNode(bigraphNodeId, offset);
		assert(params.graph.NodeOffset(nodeIndex) <= offset);
		assert(params.graph.NodeOffset(nodeIndex) + params.graph.NodeLength(nodeIndex) > offset);
		size_t offsetInNode = offset - params.graph.NodeOffset(nodeIndex);
		assert(offsetInNode < params.graph.NodeLength(nodeIndex));
		result.scores.addNodeToMap(nodeIndex);
		result.minScoreNode = nodeIndex;
		result.minScoreNodeOffset = offsetInNode;
		result.maxExactEndposScore = 0;
		result.maxExactEndposNode = result.minScoreNode;
		auto& node = result.scores.node(nodeIndex);
		node.startSlice = {0, 0, (int)offsetInNode};
		node.endSlice = {0, 0, (int)params.graph.NodeLength(nodeIndex) - 1 - (int)offsetInNode};
		node.minScore = 0;
		node.exists = true;
		for (size_t i = 1; i <= offsetInNode; i++)
		{
			size_t chunkIndex = i / (sizeof(Word) * 8);
			size_t chunkOffset = i % (sizeof(Word) * 8);
			node.HN[chunkIndex] |= ((Word)1) << chunkOffset;
		}
		for (size_t i = offsetInNode+1; i < params.graph.NodeLength(nodeIndex); i++)
		{
			size_t chunkIndex = i / (sizeof(Word) * 8);
			size_t chunkOffset = i % (sizeof(Word) * 8);
			node.HP[chunkIndex] |= ((Word)1) << chunkOffset;
		}
		return result;
	}

	static DPSlice getInitialSliceOneNodeGroup(const Params& params, const std::vector<LengthType>& nodeIndices)
	{
		DPSlice result;
		result.j = -WordConfiguration<Word>::WordSize;
		result.bandwidth = 1;
		result.minScore = 0;
		result.scores.addEmptyNodeMap(nodeIndices.size());
		for (auto nodeIndex : nodeIndices)
		{
			result.scores.addNodeToMap(nodeIndex);
			result.minScoreNode = nodeIndex;
			result.minScoreNodeOffset = params.graph.NodeLength(nodeIndex) - 1;
			auto& node = result.scores.node(nodeIndex);
			node.startSlice = {0, 0, 0};
			node.endSlice = {0, 0, 0};
			node.minScore = 0;
			node.exists = true;
		}
		return result;
	}

};

#endif
