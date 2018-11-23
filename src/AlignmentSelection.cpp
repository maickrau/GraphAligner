#include <algorithm>
#include <cstdint>
#include <cmath>
#include "AlignmentSelection.h"

//an overlap which is larger than the fraction cutoff of the smaller alignment means the alignments are incompatible
//eg alignments 12000bp and 15000bp, overlap of 12000*0.05 = 600bp means they are incompatible
const float OverlapIncompatibleFractionCutoff = 0.05;

//model the alignment as one sequence with the alphabet {match, mismatch}
//with random alignments having P(match) = P(mismatch) = 0.5
//and match having score +1 and mismatch having score -2
//then use Karlin-Altschul equation to calculate E
const float lambda = 0.48121182506; //ln(1+sqrt(5)) - ln(2)
const float K = 29.8318175495*1.25982891379;

namespace AlignmentSelection
{
	//approx ~ n(matches) - 2 * n(mismatches), higher is better
	double S(const vg::Alignment* const aln)
	{
		return aln->sequence().size() - aln->score() * 3;
	}

	double Evalue(size_t m, size_t n, const vg::Alignment* const aln)
	{
		return K * m * n * exp(-lambda * S(aln));
	}
	bool alignmentIncompatible(const vg::Alignment* const left, const vg::Alignment* const right)
	{
		auto minOverlapLen = std::min(left->sequence().size(), right->sequence().size()) * OverlapIncompatibleFractionCutoff;
		assert(left->query_position() >= 0);
		assert(right->query_position() >= 0);
		size_t leftStart = left->query_position();
		size_t leftEnd = leftStart + left->sequence().size();
		size_t rightStart = right->query_position();
		size_t rightEnd = rightStart + right->sequence().size();
		if (leftStart > rightStart)
		{
			std::swap(leftStart, rightStart);
			std::swap(leftEnd, rightEnd);
		}
		int overlap = 0;
		assert(leftStart <= rightStart);
		if (leftEnd > rightStart) overlap = leftEnd - rightStart;
		return overlap > minOverlapLen;
	}

	//lower E-value is better
	bool alignmentECompare(const vg::Alignment* const left, const vg::Alignment* const right, size_t m, size_t n)
	{
		return Evalue(m, n, left) < Evalue(m, n, right);
	}

	bool alignmentScoreCompare(const vg::Alignment* const left, const vg::Alignment* const right)
	{
		return S(left) > S(right);
	}

	//longer is better, after that lower score is better
	bool alignmentLengthCompare(const vg::Alignment* const left, const vg::Alignment* const right)
	{
		if (left->sequence().size() > right->sequence().size()) return true;
		if (right->sequence().size() > left->sequence().size()) return false;
		if (left->score() < right->score()) return true;
		return false;
	}

	std::vector<vg::Alignment> SelectAlignments(std::vector<vg::Alignment> alns, SelectionOptions options)
	{
		return SelectAlignments(alns, options, [](const vg::Alignment& aln) { return &aln; });
	}

	std::vector<vg::Alignment*> SelectAlignments(std::vector<vg::Alignment*> alns, SelectionOptions options)
	{
		return SelectAlignments(alns, options, [](vg::Alignment* aln) { return aln; });
	}
}