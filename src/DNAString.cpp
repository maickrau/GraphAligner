#include "CommonUtils.h"
#include "ThreadReadAssertion.h"
#include "DNAString.h"

std::vector<bool> getValidTwoBitChars()
{
	std::vector<bool> result;
	result.resize(256, false);
	result['a'] = true;
	result['A'] = true;
	result['c'] = true;
	result['C'] = true;
	result['g'] = true;
	result['G'] = true;
	result['t'] = true;
	result['T'] = true;
	return result;
}

std::vector<bool> getValidDNAChars()
{
	std::vector<bool> result;
	result.resize(256, false);
	result['a'] = true;
	result['A'] = true;
	result['c'] = true;
	result['C'] = true;
	result['g'] = true;
	result['G'] = true;
	result['t'] = true;
	result['T'] = true;
	result['y'] = true;
	result['Y'] = true;
	result['r'] = true;
	result['R'] = true;
	result['w'] = true;
	result['W'] = true;
	result['s'] = true;
	result['S'] = true;
	result['k'] = true;
	result['K'] = true;
	result['m'] = true;
	result['M'] = true;
	result['d'] = true;
	result['D'] = true;
	result['v'] = true;
	result['V'] = true;
	result['h'] = true;
	result['H'] = true;
	result['b'] = true;
	result['B'] = true;
	result['n'] = true;
	result['N'] = true;
	return result;
}

std::vector<size_t> getTwoBitMapping()
{
	std::vector<size_t> result;
	result.resize(256, 0);
	result['a'] = 0;
	result['A'] = 0;
	result['c'] = 1;
	result['C'] = 1;
	result['g'] = 2;
	result['G'] = 2;
	result['t'] = 3;
	result['T'] = 3;
	return result;
}

std::vector<size_t> getFourBitMapping()
{
	std::vector<size_t> result;
	result.resize(256, 0);
	result['a'] = 0;
	result['A'] = 0;
	result['c'] = 1;
	result['C'] = 1;
	result['g'] = 2;
	result['G'] = 2;
	result['t'] = 3;
	result['T'] = 3;
	result['y'] = 4;
	result['Y'] = 4;
	result['r'] = 5;
	result['R'] = 5;
	result['w'] = 6;
	result['W'] = 6;
	result['s'] = 7;
	result['S'] = 7;
	result['k'] = 8;
	result['K'] = 8;
	result['m'] = 9;
	result['M'] = 9;
	result['d'] = 10;
	result['D'] = 10;
	result['v'] = 11;
	result['V'] = 11;
	result['h'] = 12;
	result['H'] = 12;
	result['b'] = 13;
	result['B'] = 13;
	result['n'] = 14;
	result['N'] = 14;
	return result;
}

const std::vector<bool> validTwoBit = getValidTwoBitChars();
const std::vector<bool> valid = getValidDNAChars();
const std::vector<size_t> twoBitMapping = getTwoBitMapping();
const std::vector<size_t> fourBitMapping = getFourBitMapping();

DNAString::DNAString() :
storage(),
realSize(0),
lastSubstringEnd(0),
lastIndex(0),
lastOffset(0)
{
}

DNAString::DNAString(const std::string& str) :
storage(),
realSize(0),
lastSubstringEnd(0),
lastIndex(0),
lastOffset(0)
{
	storage.reserve(str.size() / 15);
	buildFromString(str);
}

DNAString& DNAString::operator=(const std::string& str)
{
	storage.clear();
	realSize = 0;
	lastSubstringEnd = 0;
	lastIndex = 0;
	lastOffset = 0;
	buildFromString(str);
	return *this;
}

void DNAString::buildFromString(const std::string& str)
{
	size_t start = 0;
	lastCharOffset = 0;
	while (start < str.size())
	{
		bool validTwoBits = true;
		for (size_t i = start; i < start+31 && i < str.size(); i++)
		{
			if (!validTwoBit[str[i]])
			{
				validTwoBits = false;
				if (!valid[str[i]])
				{
					throw CommonUtils::InvalidGraphException("Invalid sequence character: " + str[i]);
				}
				break;
			}
		}
		size_t block = 0;
		if (validTwoBits)
		{
			for (size_t i = start; i < start+31 && i < str.size(); i++)
			{
				block >>= 2;
				block += twoBitMapping[str[i]] << 60;
			}
			if (start+31 >= str.size())
			{
				lastCharOffset = str.size() - start;
				block >>= (31 - lastCharOffset) * 2;
			}
			block |= ((size_t)1 << (size_t)63);
			start += std::min((size_t)31, str.size() - start);
		}
		else
		{
			for (size_t i = start; i < start+15 && i < str.size(); i++)
			{
				block >>= 4;
				block += fourBitMapping[str[i]] << 56;
			}
			if (start+15 >= str.size())
			{
				lastCharOffset = str.size() - start;
				block >>= (15 - lastCharOffset) * 4;
			}
			start += std::min((size_t)15, str.size() - start);
		}
		storage.push_back(block);
	}
	assert(lastCharOffset > 0);
	if (storage.back() >> 63)
	{
		assert(lastCharOffset <= 31);
	}
	else
	{
		assert(lastCharOffset <= 15);
	}
	realSize = str.size();
}

size_t DNAString::size() const
{
	return realSize;
}

DNAString DNAString::reverseComplement() const
{
	assert(size() >= 1);
	assert(lastCharOffset > 0);
	size_t iterOffset = lastCharOffset - 1;
	if (storage.back() >> 63)
	{
		iterOffset *= 2;
	}
	else
	{
		iterOffset *= 4;
	}
	assert(iterOffset < 64);
	size_t iterBlock = storage.size()-1;
	DNAString result;
	std::string currentStr;
	while (true)
	{
		assert(iterOffset < 64);
		assert(iterBlock < storage.size());
		if (storage[iterBlock] >> 63)
		{
			currentStr += "TGCA"[(storage[iterBlock] >> iterOffset) & 3];
			iterOffset -= 2;
		}
		else
		{
			currentStr += "TGCARYSWMKHBDVNN"[(storage[iterBlock] >> iterOffset) & 15];
			iterOffset -= 4;
		}
		if (currentStr.size() >= 31)
		{
			size_t added = result.addString(currentStr);
			assert(added <= currentStr.size());
			currentStr.erase(currentStr.begin(), currentStr.begin() + added);
		}
		if (iterOffset > 64)
		{
			if (iterBlock == 0) break;
			iterBlock -= 1;
			if (storage[iterBlock] >> 63)
			{
				iterOffset = 60;
			}
			else
			{
				iterOffset = 56;
			}
		}
	}
	while (currentStr.size() > 0)
	{
		size_t added = result.addString(currentStr);
		assert(added <= currentStr.size());
		currentStr.erase(currentStr.begin(), currentStr.begin() + added);
	}
	assert(result.size() == size());
	return result;
}

size_t DNAString::addString(const std::string& str)
{
	bool validTwoBits = true;
	for (size_t i = 0; i < 31 && i < str.size(); i++)
	{
		if (!validTwoBit[str[i]])
		{
			validTwoBits = false;
			if (!valid[str[i]])
			{
				throw CommonUtils::InvalidGraphException("Invalid sequence character: " + str[i]);
			}
			break;
		}
	}
	size_t block = 0;
	size_t added = 0;
	if (validTwoBits)
	{
		for (size_t i = 0; i < 31 && i < str.size(); i++)
		{
			block >>= 2;
			block += twoBitMapping[str[i]] << 60;
		}
		added = std::min((size_t)31, str.size());
		if (str.size() <= 31)
		{
			lastCharOffset = str.size();
			block >>= (31 - lastCharOffset) * 2;
		}
		block |= ((size_t)1 << (size_t)63);
	}
	else
	{
		for (size_t i = 0; i < 15 && i < str.size(); i++)
		{
			block >>= 4;
			block += fourBitMapping[str[i]] << 56;
		}
		added = std::min((size_t)15, str.size());
		if (str.size() <= 15)
		{
			lastCharOffset = str.size();
			block >>= (15 - lastCharOffset) * 4;
		}
	}
	storage.push_back(block);
	realSize += added;
	return added;
}

std::string DNAString::substr(size_t start, size_t length) const
{
	assert(start == lastSubstringEnd);
	assert(start+length <= size());
	std::string result = "";
	for (size_t i = 0; i < length; i++)
	{
		if (storage[lastIndex] >> 63)
		{
			result += "ACGT"[(storage[lastIndex] >> lastOffset) & 3];
			lastOffset += 2;
			if (lastOffset == 62)
			{
				lastOffset = 0;
				lastIndex += 1;
			}
		}
		else
		{
			result += "ACGTYRWSKMDVHBNN"[(storage[lastIndex] >> lastOffset) & 15];
			lastOffset += 4;
			if (lastOffset == 60)
			{
				lastOffset = 0;
				lastIndex += 1;
			}
		}
	}
	lastSubstringEnd = start+length;
	assert(result.size() == length);
	return result;
}

std::string DNAString::toString() const
{
	lastSubstringEnd = 0;
	lastOffset = 0;
	lastIndex = 0;
	auto result = substr(0, size());
	assert(result.size() == size());
	lastSubstringEnd = 0;
	lastOffset = 0;
	lastIndex = 0;
	return result;
}

void DNAString::rewindIterators(size_t size) const
{
	assert(size > 0);
	assert(lastSubstringEnd >= size);
	for (size_t i = 0; i < size; i++)
	{
		if (storage[lastIndex] >> 63)
		{
			lastOffset -= 2;
		}
		else
		{
			lastOffset -= 4;
		}
		if (lastOffset > 64)
		{
			assert(lastIndex > 0);
			lastIndex -= 1;
			if (storage[lastIndex] >> 63)
			{
				lastOffset = 60;
			}
			else
			{
				lastOffset = 56;
			}
		}
	}
	lastSubstringEnd -= size;
}
