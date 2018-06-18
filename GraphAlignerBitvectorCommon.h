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
	using WordSlice = decltype(NodeSlice<LengthType, ScoreType, Word>::NodeSliceMapItem::startSlice);

	class EqVector
	{
	public:
		EqVector(Word BA, Word BT, Word BC, Word BG) :
		BA(BA),
		BT(BT),
		BC(BC),
		BG(BG)
		{
		}
		Word getEq(char c) const
		{
			switch(c)
			{
				case 'A':
				case 'a':
					return BA;
				case 'T':
				case 't':
					return BT;
				case 'C':
				case 'c':
					return BC;
				case 'G':
				case 'g':
					return BG;
				case '-':
				default:
					assert(false);
			}
			assert(false);
			return 0;
		}
		Word BA;
		Word BT;
		Word BC;
		Word BG;
	};

	GraphAlignerBitvectorCommon() = delete;

	static std::tuple<WordSlice, Word, Word> getNextSlice(Word Eq, WordSlice slice, Word hinP, Word hinN)
	{
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408

		Word Xv = Eq | slice.VN;
		Eq |= hinN;
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq;
		Word Ph = slice.VN | ~(Xh | slice.VP);
		Word Mh = slice.VP & Xh;
		slice.scoreEnd += Ph >> (WordConfiguration<Word>::WordSize-1);
		slice.scoreEnd -= Mh >> (WordConfiguration<Word>::WordSize-1);
		Word temp = (Mh << 1) | hinN;
		hinN = Mh >> (WordConfiguration<Word>::WordSize-1);
		slice.VP = temp | ~(Xv | Ph);
		temp = Ph << 1 | hinP;
		hinP = Ph >> (WordConfiguration<Word>::WordSize-1);
		slice.VN = temp & Xv;

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

};

#endif
