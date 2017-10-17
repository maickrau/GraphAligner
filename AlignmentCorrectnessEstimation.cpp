#include "AlignmentCorrectnessEstimation.h"
#include "ThreadReadAssertion.h"

const boost::rational<boost::multiprecision::cpp_int> AlignmentCorrectnessEstimationState::correctMismatchProbability {15, 100}; //15% from pacbio error rate
const boost::rational<boost::multiprecision::cpp_int> AlignmentCorrectnessEstimationState::falseMismatchProbability {50, 100}; //50% empirically
const boost::rational<boost::multiprecision::cpp_int> AlignmentCorrectnessEstimationState::falseToCorrectTransitionProbability {1, 100000}; //10^-5. one slice at <=12 mismatches or two slices at <=15 mismatchs
const boost::rational<boost::multiprecision::cpp_int> AlignmentCorrectnessEstimationState::correctToFalseTransitionProbability {1LL, 1000000000000000LL}; //10^-15. three slices at >=25 mismatches or two slices at >=30 mismatches

AlignmentCorrectnessEstimationState::AlignmentCorrectnessEstimationState() :
correctProbability(80, 100), //80% arbitrarily
falseProbability(20, 100), //20% arbitrarily
correctFromCorrectTrace(false),
falseFromCorrectTrace(false)
{
}

bool AlignmentCorrectnessEstimationState::CurrentlyCorrect() const
{
	return correctProbability > falseProbability;
}

bool AlignmentCorrectnessEstimationState::CorrectFromCorrect() const
{
	return correctFromCorrectTrace;
}

bool AlignmentCorrectnessEstimationState::FalseFromCorrect() const
{
	return falseFromCorrectTrace;
}

template <typename T>
T factorial(T n)
{
	T result {1};
	for (int i = 2; i <= n; i += 1)
	{
		result *= i;
	}
	return result;
}

template <typename T>
T choose(T n, T k)
{
	return factorial<T>(n) / factorial<T>(k) / factorial<T>(n - k);
}

template <typename T>
T powr(T base, int exponent)
{
	assert(exponent >= 0);
	if (exponent == 0) return T{1};
	if (exponent == 1) return base;
	if (exponent % 2 == 0)
	{
		auto part = powr(base, exponent / 2);
		assert(part > 0);
		return part * part;
	}
	if (exponent % 2 == 1)
	{
		auto part = powr(base, exponent / 2);
		assert(part > 0);
		return part * part * base;
	}
	assert(false);
	std::abort();
	return T{};
}

AlignmentCorrectnessEstimationState AlignmentCorrectnessEstimationState::NextState(int mismatches, int rowSize) const
{
	assert(mismatches >= 0);
	assert(mismatches <= rowSize);
	AlignmentCorrectnessEstimationState result;
	result.correctFromCorrectTrace = correctProbability * (1 - correctToFalseTransitionProbability) >= falseProbability * falseToCorrectTransitionProbability;
	result.falseFromCorrectTrace = correctProbability * correctToFalseTransitionProbability >= falseProbability * (1 - falseToCorrectTransitionProbability);
	boost::rational<boost::multiprecision::cpp_int> newCorrectProbability = std::max(correctProbability * (1 - correctToFalseTransitionProbability), falseProbability * falseToCorrectTransitionProbability);
	boost::rational<boost::multiprecision::cpp_int> newFalseProbability = std::max(correctProbability * correctToFalseTransitionProbability, falseProbability * (1 - falseToCorrectTransitionProbability));
	auto chooseresult = choose<boost::multiprecision::cpp_int>(rowSize, mismatches);
	auto correctMultiplier = chooseresult * powr(correctMismatchProbability, mismatches) * powr(1 - correctMismatchProbability, rowSize - mismatches);
	auto falseMultiplier = chooseresult * powr(falseMismatchProbability, mismatches) * powr(1 - falseMismatchProbability, rowSize - mismatches);
	newCorrectProbability *= correctMultiplier;
	newFalseProbability *= falseMultiplier;
	result.correctProbability = newCorrectProbability;
	result.falseProbability = newFalseProbability;
	auto normalizer = result.correctProbability + result.falseProbability;
	result.correctProbability /= normalizer;
	result.falseProbability /= normalizer;
	return result;
}
