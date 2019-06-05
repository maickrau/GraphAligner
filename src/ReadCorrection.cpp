#include "ThreadReadAssertion.h"
#include "ReadCorrection.h"

std::string toUpper(std::string seq)
{
	for (auto& c : seq)
	{
		c = toupper(c);
	}
	return seq;
}

std::string toLower(std::string seq)
{
	for (auto& c : seq)
	{
		c = tolower(c);
	}
	return seq;
}

size_t getLongestOverlap(const std::string& left, const std::string& right, size_t maxOverlap)
{
	if (left.size() < maxOverlap) maxOverlap = left.size();
	if (right.size() < maxOverlap) maxOverlap = right.size();
	for (size_t i = maxOverlap; i > 0; i--)
	{
		bool match = true;
		for (size_t a = 0; a < i && match; a++)
		{
			if (left[left.size() - maxOverlap + a] != right[a]) match = false;
		}
		if (match) return i;
	}
	return 0;
}

std::string getCorrected(const std::string& raw, const std::vector<Correction>& corrections, size_t maxOverlap)
{
	std::string result;
	size_t currentEnd = 0;
	for (size_t i = 0; i < corrections.size(); i++)
	{
		assert(i == 0 || corrections[i].startIndex >= corrections[i-1].startIndex);
		if (corrections[i].startIndex < currentEnd)
		{
			size_t overlap = getLongestOverlap(result, corrections[i].corrected, maxOverlap);
			result += toUpper(corrections[i].corrected.substr(overlap));
		}
		else if (corrections[i].startIndex > currentEnd)
		{
			result += toLower(raw.substr(currentEnd, corrections[i].startIndex - currentEnd));
			result += toUpper(corrections[i].corrected);
		}
		else
		{
			assert(corrections[i].startIndex == currentEnd);
			result += toUpper(corrections[i].corrected);
		}
		currentEnd = corrections[i].endIndex;
	}
	if (currentEnd < raw.size()) result += toLower(raw.substr(currentEnd));
	return result;
}
