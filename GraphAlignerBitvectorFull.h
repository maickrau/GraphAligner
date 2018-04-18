#ifndef GraphAlignerBitvectorFull_h
#define GraphAlignerBitvectorFull_h

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerBitvectorFull
{
private:
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using BV = GraphAlignerBitvectorCommon<LengthType, ScoreType, Word>;
	using EqVector = typename BV::EqVector;
	using Params = typename Common::Params;
	const Params params;
public:
	GraphAlignerBitvectorFull(const Params& params) :
	params(params)
	{
	}
	size_t alignAndGetScore(const std::string& sequence) const
	{
		std::vector<WordSlice> currentSlice;
		std::vector<WordSlice> previousSlice;
		currentSlice.resize(params.graph.NodeSequencesSize() { WordConfiguration<Word>::AllOnes, WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::WordSize, 0, true });
		previousSlice.resize(params.graph.NodeSequencesSize() { WordConfiguration<Word>::AllZeros, WordConfiguration<Word>::AllZeros, 0, 0, true });
		ArrayPriorityQueue<NodeWithPriority> calculableQueue { sequence.size() };
		int slices = (sequence.size() + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		for (int slice = 0; slice < slices.size(); slices++)
		{
			size_t j = slice * WordConfiguration<Word>::WordSize;
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

			while (calculableQueue.size() > 0)
			{
				
			}
			std::swap(previousSlice, currentSlice);
		}
	}
}

#endif
