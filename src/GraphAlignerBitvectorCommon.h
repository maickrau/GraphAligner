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
#include "AlignmentCorrectnessEstimation.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"
#include "GraphAlignerCommon.h"
#include "ArrayPriorityQueue.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorCommon
{
private:
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using Params = typename Common::Params;
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

	GraphAlignerBitvectorCommon() = delete;

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static std::tuple<WordSlice, Word, Word> getNextSlice(Word Eq, WordSlice slice, Word hinP, Word hinN)
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

		return std::make_tuple(slice, hinP, hinN);
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
		Word BA = WordConfiguration<Word>::AllZeros;
		Word BT = WordConfiguration<Word>::AllZeros;
		Word BC = WordConfiguration<Word>::AllZeros;
		Word BG = WordConfiguration<Word>::AllZeros;
		for (int i = 0; i < WordConfiguration<Word>::WordSize && j+i < sequence.size(); i++)
		{
			Word mask = ((Word)1) << i;
			if (Common::characterMatch(sequence[j+i], 'A')) BA |= mask;
			if (Common::characterMatch(sequence[j+i], 'C')) BC |= mask;
			if (Common::characterMatch(sequence[j+i], 'T')) BT |= mask;
			if (Common::characterMatch(sequence[j+i], 'G')) BG |= mask;
		}
		if (j + WordConfiguration<Word>::WordSize > sequence.size())
		{
			Word mask = WordConfiguration<Word>::AllOnes << (WordConfiguration<Word>::WordSize - j + sequence.size());
			assert((BA & mask) == 0);
			assert((BT & mask) == 0);
			assert((BC & mask) == 0);
			assert((BG & mask) == 0);
			BA |= mask;
			BT |= mask;
			BC |= mask;
			BG |= mask;
		}
		assert((BA | BC | BT | BG) == WordConfiguration<Word>::AllOnes);
		EqVector EqV {BA, BT, BC, BG};
		return EqV;
	}

};

#endif
