#ifndef AlignmentCorrectnessEstimation_h
#define AlignmentCorrectnessEstimation_h

#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>

class AlignmentCorrectnessEstimationState
{
public:
	AlignmentCorrectnessEstimationState();
	bool CurrentlyCorrect() const;
	bool CorrectFromCorrect() const;
	bool FalseFromCorrect() const;
	AlignmentCorrectnessEstimationState NextState(int mismatches, int rowSize) const;
private:
	const static boost::rational<boost::multiprecision::cpp_int> correctMismatchProbability;
	const static boost::rational<boost::multiprecision::cpp_int> falseMismatchProbability;
	const static boost::rational<boost::multiprecision::cpp_int> falseToCorrectTransitionProbability;
	const static boost::rational<boost::multiprecision::cpp_int> correctToFalseTransitionProbability;
	boost::rational<boost::multiprecision::cpp_int> correctProbability;
	boost::rational<boost::multiprecision::cpp_int> falseProbability;
	bool correctFromCorrectTrace;
	bool falseFromCorrectTrace;
};

#endif
