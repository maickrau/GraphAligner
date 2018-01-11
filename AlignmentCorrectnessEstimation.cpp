#include <cmath>
#include <vector>
#include "AlignmentCorrectnessEstimation.h"
#include "ThreadReadAssertion.h"

const double correctMismatchLogProbability = log(0.15); //15% from pacbio error rate
const double correctMatchLogProbability = log(1.0 - 0.15);
const double falseMismatchLogProbability = log(0.5); //50% empirically
const double falseMatchLogProbability = log(1.0 - 0.5);
const double falseToCorrectTransitionLogProbability = log(0.00001); //10^-5. one slice at <=12 mismatches or two slices at <=15 mismatchs
const double falseToFalseTransitionLogProbability = log(1.0 - 0.00001);
const double correctToFalseTransitionLogProbability = log(0.000000000000001); //10^-15. three slices at >=25 mismatches or two slices at >=30 mismatches
const double correctToCorrectTransitionLogProbability = log(1.0 - 0.000000000000001);

std::vector<double> getLogFactorials()
{
	std::vector<double> result;
	result.push_back(0);
	for (int i = 1; i <= 64; i++)
	{
		result.push_back(result.back() + log(i));
	}
	return result;
}

const std::vector<double> logFactorials = getLogFactorials();

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

double logChoose(int n, int k)
{
	return logFactorials[n] - logFactorials[k] - logFactorials[n - k];
}

double logPowr(double logBase, int exponent)
{
	return exponent * logBase;
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
	assert(rowSize == 64 || rowSize == 1);
	assert(mismatches >= 0);
	assert(mismatches <= rowSize);
	AlignmentCorrectnessEstimationState result;
	result.correctFromCorrectTrace = correctLogOdds + correctToCorrectTransitionLogProbability >= falseLogOdds + falseToCorrectTransitionLogProbability;
	result.falseFromCorrectTrace = correctLogOdds + correctToFalseTransitionLogProbability >= falseLogOdds + falseToFalseTransitionLogProbability;
	double newCorrectProbability = std::max(correctLogOdds + correctToCorrectTransitionLogProbability, falseLogOdds + falseToCorrectTransitionLogProbability);
	double newFalseProbability = std::max(correctLogOdds + correctToFalseTransitionLogProbability, falseLogOdds + falseToFalseTransitionLogProbability);
	auto chooseresult = logChoose(rowSize, mismatches);
	auto correctMultiplier = chooseresult + logPowr(correctMismatchLogProbability, mismatches) + logPowr(correctMatchLogProbability, rowSize - mismatches);
	auto falseMultiplier = chooseresult + logPowr(falseMismatchLogProbability, mismatches) + logPowr(falseMatchLogProbability, rowSize - mismatches);
	newCorrectProbability += correctMultiplier;
	newFalseProbability += falseMultiplier;
	result.correctLogOdds = newCorrectProbability;
	result.falseLogOdds = newFalseProbability;
	return result;
}
