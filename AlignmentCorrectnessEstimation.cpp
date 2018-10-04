#include <cmath>
#include <vector>
#include "AlignmentCorrectnessEstimation.h"
#include "ThreadReadAssertion.h"

//empirically from aligning one ONT to its correct position in the genome
const double correctMean = 0.1875;
const double correctStddev = 0.0955;
//empirically from aligning one random read to a position
const double wrongMean = 0.5;
const double wrondStddev = 0.0291;

const int wordSize = 64;

const double falseToCorrectTransitionLogProbability = log(0.00001); //10^-5. arbitrary.
const double falseToFalseTransitionLogProbability = log(1.0 - 0.00001);
const double correctToFalseTransitionLogProbability = log(0.0000000001); //10^-10. arbitrary.
const double correctToCorrectTransitionLogProbability = log(1.0 - 0.0000000001);

double stddistlog(double val, double mean, double stddev)
{
	return -(val-mean)*(val-mean)/(2*stddev*stddev);
}

void normalize(std::vector<double>& logs)
{
	double sum = 0;
	for (auto x : logs)
	{
		sum += exp(x);
	}
	double add = log(1.0/sum);
	for (auto& x : logs)
	{
		x += add;
	}
}

std::vector<double> getCorrectLogOdds()
{
	std::vector<double> result;
	for (int i = 0; i <= wordSize/2; i++)
	{
		result.push_back(stddistlog(i, correctMean*wordSize, correctStddev*wordSize));
	}
	normalize(result);
	for (int i = wordSize/2; i < wordSize; i++)
	{
		result.push_back(result.back());
	}
	return result;
}

std::vector<double> getWrongLogOdds()
{
	std::vector<double> result;
	for (int i = 0; i <= wordSize/2; i++)
	{
		result.push_back(stddistlog(i, wrongMean*wordSize, wrondStddev*wordSize));
	}
	normalize(result);
	for (int i = wordSize/2; i < wordSize; i++)
	{
		result.push_back(result.back());
	}
	return result;
}

const std::vector<double> precomputedCorrectLogOdds = getCorrectLogOdds();
const std::vector<double> precomputedWrongLogOdds = getWrongLogOdds();

AlignmentCorrectnessEstimationState::AlignmentCorrectnessEstimationState() :
correctLogOdds(log(0.8)), //80% arbitrarily
falseLogOdds(log(0.2)), //20% arbitrarily
correctFromCorrectTrace(false),
falseFromCorrectTrace(false)
{
}

bool AlignmentCorrectnessEstimationState::CurrentlyCorrect() const
{
	return correctLogOdds > falseLogOdds;
}

bool AlignmentCorrectnessEstimationState::CorrectFromCorrect() const
{
	return correctFromCorrectTrace;
}

bool AlignmentCorrectnessEstimationState::FalseFromCorrect() const
{
	return falseFromCorrectTrace;
}

double AlignmentCorrectnessEstimationState::CorrectLogOdds() const
{
	return correctLogOdds;
}

double AlignmentCorrectnessEstimationState::FalseLogOdds() const
{
	return falseLogOdds;
}

AlignmentCorrectnessEstimationState AlignmentCorrectnessEstimationState::NextState(int mismatches, int rowSize) const
{
	assert(rowSize == 64);
	// assert(rowSize == 64 || rowSize == 1);
	assert(mismatches >= 0);
	AlignmentCorrectnessEstimationState result;
	result.correctFromCorrectTrace = correctLogOdds + correctToCorrectTransitionLogProbability >= falseLogOdds + falseToCorrectTransitionLogProbability;
	result.falseFromCorrectTrace = correctLogOdds + correctToFalseTransitionLogProbability >= falseLogOdds + falseToFalseTransitionLogProbability;
	double newCorrectProbability = std::max(correctLogOdds + correctToCorrectTransitionLogProbability, falseLogOdds + falseToCorrectTransitionLogProbability);
	double newFalseProbability = std::max(correctLogOdds + correctToFalseTransitionLogProbability, falseLogOdds + falseToFalseTransitionLogProbability);
	assert(precomputedCorrectLogOdds.size() == precomputedWrongLogOdds.size());
	if (mismatches < precomputedCorrectLogOdds.size())
	{
		newCorrectProbability += precomputedCorrectLogOdds[mismatches];
		newFalseProbability += precomputedWrongLogOdds[mismatches];
	}
	else
	{
		newCorrectProbability += precomputedCorrectLogOdds.back();
		newFalseProbability += precomputedWrongLogOdds.back();
	}
	result.correctLogOdds = newCorrectProbability;
	result.falseLogOdds = newFalseProbability;
	return result;
}
