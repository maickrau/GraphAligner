#include <vector>
#include <tuple>
#include "ThreadReadAssertion.h"
#include "ByteStuff.h"

size_t index(size_t scorediffPlusSeventeen, size_t sign, size_t low, size_t high)
{
	assert(scorediffPlusSeventeen <= 17+17);
	assert(sign < 256);
	assert(low < 256);
	assert(high < 256);
	return scorediffPlusSeventeen*256*256*256+sign*256*256+low*256+high;
}

std::tuple<uint8_t, uint8_t, int8_t> precalcVPVNChange(int8_t scorediff, uint8_t sign, uint8_t low, uint8_t high)
{
	uint8_t leftSmaller = 0;
	uint8_t rightSmaller = 0;
	int8_t leftscore = 0;
	int8_t rightscore = scorediff;
	for (int i = 0; i < 8; i++)
	{
		uint8_t mask = 1 << i;
		if (sign & mask)
		{
			if (low & mask) leftscore--;
			if (high & mask) leftscore--;
		}
		else
		{
			if (low & mask) rightscore--;
			if (high & mask) rightscore--;
		}
		if (leftscore < rightscore) leftSmaller |= mask;
		if (rightscore < leftscore) rightSmaller |= mask;
	}
	return std::make_tuple(leftSmaller, rightSmaller, rightscore - leftscore - scorediff);
}

std::vector<std::tuple<uint8_t, uint8_t, int8_t>> getPrecalcedChanges()
{
	std::vector<std::tuple<uint8_t, uint8_t, int8_t>> result;
	result.reserve((17+17+1)*256*256*256);
	for (int scorediff = -17; scorediff <= 17; scorediff++)
	{
		for (int sign = 0; sign < 256; sign++)
		{
			for (int low = 0; low < 256; low++)
			{
				for (int high = 0; high < 256; high++)
				{
					assert(result.size() == index(scorediff, sign, low, high));
					result.push_back(precalcVPVNChange(scorediff, sign, low, high));
				}
			}
		}
	}
	assert(result.size() == (17+17+1)*256*256*256);
	return result;
}

std::vector<std::tuple<uint8_t, uint8_t, int8_t>> precalcedVPVNChanges = getPrecalcedChanges();

namespace ByteStuff
{
	std::tuple<uint8_t, uint8_t, int8_t> VPVNChange(size_t scorediffPlusSeventeen, size_t sign, size_t low, size_t high)
	{
		size_t vecindex = index(scorediffPlusSeventeen, sign, low, high);
		assert(vecindex < precalcedVPVNChanges.size());
		return precalcedVPVNChanges[vecindex];
	}
}
