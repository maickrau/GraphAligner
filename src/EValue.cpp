#include <vector>
#include <cmath>
#include "ThreadReadAssertion.h"
#include "EValue.h"

// model the alignment as one sequence with the alphabet {match, mismatch}
// with random alignments having P(match) = P(mismatch) = 0.5
// and match having score +1, mismatch score <0 chosen such that an alignment at minIdentity has score 0
// then use Karlin-Altschul equation to calculate E
// get lambda and K numerically depending no the match & mismatch score

// ...except the bitvector algorithm doesn't give the number of matches and mismatches
// so approximate the alignment score as (length * matchScore + numEdits * (mismatchScore-matchScore))
// this is close enough hopefully

constexpr double e = 2.71828182845904523536028747135266249775724709369995;

EValueCalculator::EValueCalculator() :
matchScore(-1),
mismatchScore(-1),
lambda(-1),
K(-1)
{
}

EValueCalculator::EValueCalculator(double minIdentity) :
matchScore(1),
mismatchScore(-minIdentity / (1.0 - minIdentity)),
lambda(-1),
K(-1)
{
	initializeLambda();
	initializeK();
}

double EValueCalculator::getEValue(size_t databaseSize, size_t querySize, double alignmentScore) const
{
	return K * databaseSize * querySize * pow(e, -lambda * alignmentScore);
}

double EValueCalculator::getEValue(size_t databaseSize, size_t querySize, size_t alignmentLength, size_t numEdits) const
{
	return getEValue(databaseSize, querySize, getAlignmentScore(alignmentLength, numEdits));
}

double EValueCalculator::getAlignmentScore(size_t alignmentLength, size_t numEdits) const
{
	return alignmentLength * matchScore - numEdits * (mismatchScore - matchScore);
}

void EValueCalculator::initializeLambda()
{
	// lambda is bounded by 0 < lambda < ln(2) < 0.7
	double guessMin = 0;
	double guessMax = 0.7;
	// bisect, max error 2^-100
	for (int i = 0; i < 100; i++)
	{
		double guessMid = (guessMin + guessMax) * 0.5;
		double valueMid = pow(e, guessMid*matchScore) * .5 + pow(e, guessMid*mismatchScore) * 0.5 - 1;
		if (valueMid < 0) guessMin = guessMid;
		if (valueMid > 0) guessMax = guessMid;
		// due to floating point precision limits
		if (valueMid == 0)
		{
			guessMin = guessMid;
			guessMax = guessMid;
			break;
		}
		if (guessMin == guessMax) break;
	}
	lambda = (guessMin + guessMax) / 2;
}

void EValueCalculator::initializeK()
{
	assert(lambda != -1);
	double seriesSum = 0;
	std::vector<size_t> pascalsTriangle;
	pascalsTriangle.push_back(1);
	for (int k = 1; k < 10; k++)
	{
		std::vector<size_t> newTriangle;
		newTriangle.resize(pascalsTriangle.size()+1, 0);
		for (size_t j = 0; j < pascalsTriangle.size(); j++)
		{
			newTriangle[j] += pascalsTriangle[j];
			newTriangle[j+1] += pascalsTriangle[j];
		}
		pascalsTriangle = newTriangle;
		assert(pascalsTriangle[0] == 1);
		assert(pascalsTriangle.back() == 1);
		assert(pascalsTriangle.size() == (size_t)k+1);
		size_t triangleSum = 0;
		for (auto n : pascalsTriangle) triangleSum += n;
		double negativeExpectation = 0;
		double greaterProbability = 0;
		for (size_t j = 0; j < pascalsTriangle.size(); j++)
		{
			size_t matches = j;
			size_t mismatches = pascalsTriangle.size() - 1 - j;
			double score = (double)matches * matchScore + (double)mismatches * mismatchScore;
			double probability = (double)pascalsTriangle[j] / (double)triangleSum;
			if (score < 0) negativeExpectation += pow(e, lambda * score) * probability;
			if (score >= 0) greaterProbability += probability;
		}
		seriesSum += (negativeExpectation + greaterProbability) / (double)k;
	}
	double expectation = .5 * matchScore * pow(e, lambda * matchScore) + .5 * mismatchScore * pow(e, lambda * mismatchScore);
	double Cstar = pow(e, -2 * seriesSum) / (lambda * expectation);
	// assume delta is 1 even though its not really
	K = Cstar * lambda / (1.0 - pow(e, -lambda));
}
