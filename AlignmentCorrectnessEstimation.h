#ifndef AlignmentCorrectnessEstimation_h
#define AlignmentCorrectnessEstimation_h

class AlignmentCorrectnessEstimationState
{
public:
	AlignmentCorrectnessEstimationState();
	bool CurrentlyCorrect() const;
	bool CorrectFromCorrect() const;
	bool FalseFromCorrect() const;
	double CorrectLogOdds() const;
	double FalseLogOdds() const;
	AlignmentCorrectnessEstimationState NextState(int mismatches, int rowSize) const;
private:
	double correctLogOdds;
	double falseLogOdds;
	bool correctFromCorrectTrace;
	bool falseFromCorrectTrace;
};

#endif
