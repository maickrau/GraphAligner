#ifndef EValue_h
#define EValue_h

class EValueCalculator
{
public:
	EValueCalculator();
	EValueCalculator(double minIdentity);
	double getAlignmentScore(size_t alignmentLength, size_t numEdits) const;
	double getEValue(size_t databaseSize, size_t querySize, size_t alignmentLength, size_t numEdits) const;
	double getEValue(size_t databaseSize, size_t querySize, double alignmentScore) const;
private:
	void initializeLambda();
	void initializeK();
	double matchScore;
	double mismatchScore;
	double lambda;
	double K;
};

#endif