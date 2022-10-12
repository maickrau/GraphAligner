#include <cassert>
#include "RankBitvector.h"
#include "Serialize.h"

int popcount(uint64_t x)
{
	//https://gcc.gnu.org/onlinedocs/gcc-4.8.4/gcc/X86-Built-in-Functions.html
	// return __builtin_popcountll(x);
	//for some reason __builtin_popcount takes 21 instructions so call assembly directly
	__asm__("popcnt %0, %0" : "+r" (x));
	return x;
}

RankBitvector::RankBitvector() :
values(),
ranksBuilt(false),
realSize(0)
{
}

RankBitvector::RankBitvector(size_t size) :
values(),
ranksBuilt(false),
realSize(0)
{
	resize(size);
}

void RankBitvector::resize(size_t newSize)
{
	assert(!ranksBuilt);
	realSize = newSize;
	values.resize(10 * ((newSize + 511) / 512), 0);
}

void RankBitvector::clear()
{
	assert(!ranksBuilt);
	for (size_t i = 0; i < values.size(); i++)
	{
		values[i] = 0;
	}
}

void RankBitvector::buildRanks()
{
	assert(!ranksBuilt);
	uint64_t bigBlockSum = 0;
	for (size_t start = 0; start < values.size(); start += 10)
	{
		assert(values[start] == 0);
		assert(values[start+1] == 0);
		values[start] = bigBlockSum;
		uint64_t smallBlockSum = 0;
		for (size_t i = 0; i < 7; i++)
		{
			smallBlockSum += popcount(values[start+2+i]);
			assert(smallBlockSum < 512);
			values[start+1] |= smallBlockSum << (uint64_t)(i * 9);
		}
		smallBlockSum += popcount(values[start+2+7]);
		bigBlockSum += smallBlockSum;
	}
	ranksBuilt = true;
}

size_t RankBitvector::size() const
{
	return realSize;
}

bool RankBitvector::get(size_t index) const
{
	assert(index < size());
	size_t bigBlock = index / 512;
	size_t blockStart = bigBlock * 10;
	size_t smallBlockIndex = blockStart + 2 + (index % 512) / 64;
	size_t offset = index % 64;
	uint64_t value = values[smallBlockIndex];
	return (bool)(((value >> offset) & 1) == 1);
}

void RankBitvector::set(size_t index, bool value)
{
	assert(index < size());
	assert(!ranksBuilt);
	size_t bigBlock = index / 512;
	size_t blockStart = bigBlock * 10;
	size_t smallBlockIndex = blockStart + 2 + (index % 512) / 64;
	size_t offset = index % 64;
	uint64_t oldValue = values[smallBlockIndex];
	if (value)
	{
		oldValue |= (uint64_t)1 << (uint64_t)offset;
	}
	else
	{
		oldValue &= ~((uint64_t)1 << (uint64_t)offset);
	}
	values[smallBlockIndex] = oldValue;
}

size_t RankBitvector::rankZero(size_t index) const
{
	return index - rankOne(index);
}

size_t RankBitvector::rankOne(size_t index) const
{
	assert(index <= size());
	assert(ranksBuilt);
	size_t bigBlock = index / 512;
	size_t blockStart = bigBlock * 10;
	size_t smallBlock = (index % 512) / 64;
	assert(smallBlock <= 8);
	size_t smallBlockIndex = blockStart + 2 + smallBlock;
	size_t offset = index % 64;
	size_t result = values[blockStart];
	if (smallBlock > 0)
	{
		uint64_t smallResult = (values[blockStart+1] >> (uint64_t)((smallBlock-1) * 9)) & (uint64_t)511;
		result += smallResult;
	}
	uint64_t bitBlock = values[smallBlockIndex];
	result += popcount(bitBlock & (((uint64_t)1 << (uint64_t)(offset)) - 1));
	return result;
}

void RankBitvector::save(std::ostream& stream) const
{
	assert(ranksBuilt);
	serialize(stream, values);
	serialize(stream, realSize);
}

void RankBitvector::load(std::istream& stream)
{
	assert(!ranksBuilt);
	deserialize(stream, values);
	deserialize(stream, realSize);
	ranksBuilt = true;
}

bool RankBitvector::operator==(const RankBitvector& other) const
{
	if (realSize != other.realSize) return false;
	if (ranksBuilt != other.ranksBuilt) return false;
	if (values != other.values) return false;
	return true;
}

bool RankBitvector::operator!=(const RankBitvector& other) const
{
	return !(*this == other);
}

bool RankBitvector::operator[](size_t index) const
{
	return get(index);
}
