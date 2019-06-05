#ifndef ReadCorrection_h
#define ReadCorrection_h

#include <string>
#include <vector>

struct Correction
{
	size_t startIndex;
	size_t endIndex;
	std::string corrected;
};

std::string getCorrected(const std::string& raw, const std::vector<Correction>& corrections, size_t maxOverlap);

#endif
