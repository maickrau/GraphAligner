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
	using WordSlice = typename WordContainer<LengthType, ScoreType, Word>::Slice;
public:
	class NodeWithPriority
	{
	public:
		NodeWithPriority(LengthType node, size_t offset, size_t endOffset, int priority) : node(node), offset(offset), endOffset(endOffset), priority(priority) {}
		bool operator>(const NodeWithPriority& other) const
		{
		return priority > other.priority;
		}
		bool operator<(const NodeWithPriority& other) const
		{
		return priority < other.priority;
		}
		LengthType node;
		size_t offset;
		size_t endOffset;
		int priority;
	};

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

	static EqVector getEqVector(const std::string& sequence, size_t pos)
	{
		//preprocessed bitvectors for character equality
		Word BA = WordConfiguration<Word>::AllZeros;
		Word BT = WordConfiguration<Word>::AllZeros;
		Word BC = WordConfiguration<Word>::AllZeros;
		Word BG = WordConfiguration<Word>::AllZeros;
		for (int i = 0; i < WordConfiguration<Word>::WordSize && pos+i < sequence.size(); i++)
		{
			Word mask = ((Word)1) << i;
			if (Common::characterMatch(sequence[pos+i], 'A')) BA |= mask;
			if (Common::characterMatch(sequence[pos+i], 'C')) BC |= mask;
			if (Common::characterMatch(sequence[pos+i], 'T')) BT |= mask;
			if (Common::characterMatch(sequence[pos+i], 'G')) BG |= mask;
		}
		if (pos + WordConfiguration<Word>::WordSize > sequence.size())
		{
			Word mask = WordConfiguration<Word>::AllOnes << (WordConfiguration<Word>::WordSize - pos + sequence.size());
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
		return EqVector {BA, BT, BC, BG};
	}

	static void assertSliceCorrectness(const WordSlice& current, const WordSlice& up, bool previousBand)
	{
	#ifndef NDEBUG
		auto wcvp = WordConfiguration<Word>::popcount(current.VP);
		auto wcvn = WordConfiguration<Word>::popcount(current.VN);
		assert(current.scoreEnd == current.scoreBeforeStart + wcvp - wcvn);

		assert(current.scoreBeforeStart >= 0);
		assert(current.scoreEnd >= 0);
		assert(current.scoreBeforeStart <= current.scoreEnd + WordConfiguration<Word>::WordSize);
		assert(current.scoreEnd <= current.scoreBeforeStart + WordConfiguration<Word>::WordSize);
		assert((current.VP & current.VN) == WordConfiguration<Word>::AllZeros);

		assert(!previousBand || current.scoreBeforeStart <= up.scoreEnd);
		assert(current.scoreBeforeStart >= 0);
	#endif
	}

	GraphAlignerBitvectorCommon() = delete;

	static WordSlice getNextSlice(Word Eq, WordSlice slice, bool upInsideBand, bool upleftInsideBand, bool diagonalInsideBand, bool previousEq, WordSlice previous, WordSlice up)
	{
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408

		auto oldValue = slice.scoreBeforeStart;
		if (!slice.scoreBeforeExists) Eq &= ~((Word)1);
		slice.scoreBeforeExists = upInsideBand;
		if (!diagonalInsideBand) Eq &= ~((Word)1);
		if (!upleftInsideBand)
		{
			slice.scoreBeforeStart += 1;
		}
		else
		{
			const auto lastBitMask = ((Word)1) << (WordConfiguration<Word>::WordSize - 1);
			assert(slice.scoreBeforeStart <= previous.scoreEnd);
			slice.scoreBeforeStart = std::min(slice.scoreBeforeStart + 1, previous.scoreEnd - ((previous.VP & lastBitMask) ? 1 : 0) + ((previous.VN & lastBitMask) ? 1 : 0) + (previousEq ? 0 : 1));
		}
		auto hin = slice.scoreBeforeStart - oldValue;

		Word Xv = Eq | slice.VN;
		//between 7 and 8
		if (hin < 0) Eq |= 1;
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq;
		Word Ph = slice.VN | ~(Xh | slice.VP);
		Word Mh = slice.VP & Xh;
		// if (~Ph & confirmedMask) confirmTentative = true;
		const Word lastBitMask = (((Word)1) << (WordConfiguration<Word>::WordSize - 1));
		if (Ph & lastBitMask)
		{
			slice.scoreEnd += 1;
		}
		else if (Mh & lastBitMask)
		{
			slice.scoreEnd -= 1;
		}
		Ph <<= 1;
		Mh <<= 1;
		//between 16 and 17
		if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
		slice.VP = Mh | ~(Xv | Ph);
		slice.VN = Ph & Xv;

		slice.sliceExists = true;

#ifndef NDEBUG
		auto wcvp = WordConfiguration<Word>::popcount(slice.VP);
		auto wcvn = WordConfiguration<Word>::popcount(slice.VN);
		assert(slice.scoreEnd == slice.scoreBeforeStart + wcvp - wcvn);
#endif

		return slice;
	}

	static WordSlice getNextSliceFullBand(Word Eq, WordSlice slice, bool previousEq, WordSlice previous)
	{
		//http://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf
		//pages 405 and 408

		auto oldValue = slice.scoreBeforeStart;
		const Word lastBitMask = ((Word)1) << (WordConfiguration<Word>::WordSize - 1);
		assert(slice.scoreBeforeStart <= previous.scoreEnd);
		slice.scoreBeforeStart = std::min(slice.scoreBeforeStart + 1, previous.scoreEnd - ((previous.VP & lastBitMask) ? 1 : 0) + ((previous.VN & lastBitMask) ? 1 : 0) + (previousEq ? 0 : 1));
		auto hin = slice.scoreBeforeStart - oldValue;

		Word Xv = Eq | slice.VN;
		//between 7 and 8
		if (hin < 0) Eq |= 1;
		Word Xh = (((Eq & slice.VP) + slice.VP) ^ slice.VP) | Eq;
		Word Ph = slice.VN | ~(Xh | slice.VP);
		Word Mh = slice.VP & Xh;
		if (Ph & lastBitMask)
		{
			slice.scoreEnd += 1;
		}
		else if (Mh & lastBitMask)
		{
			slice.scoreEnd -= 1;
		}
		Ph <<= 1;
		Mh <<= 1;
		//between 16 and 17
		if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
		slice.VP = Mh | ~(Xv | Ph);
		slice.VN = Ph & Xv;

#ifndef NDEBUG
		auto wcvp = WordConfiguration<Word>::popcount(slice.VP);
		auto wcvn = WordConfiguration<Word>::popcount(slice.VN);
		assert(slice.scoreEnd == slice.scoreBeforeStart + wcvp - wcvn);
#endif

		return slice;
	}

	static WordSlice flattenWordSlice(WordSlice slice, size_t row)
	{
		Word mask = ~(WordConfiguration<Word>::AllOnes << row);
		slice.VP &= mask;
		slice.VN &= mask;
		slice.scoreEnd = slice.scoreBeforeStart + WordConfiguration<Word>::popcount(slice.VP) - WordConfiguration<Word>::popcount(slice.VN);
		return slice;
	}

};

#endif
