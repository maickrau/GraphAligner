#ifndef DNAString_h
#define DNAString_h

#include <vector>
#include "ThreadReadAssertion.h"

template <typename Word>
class DNAString
{
	static constexpr size_t BITSIZE = sizeof(Word)*8;
	static constexpr size_t CHARSINBLOCK = BITSIZE / 2;
public:
	DNAString() : 
	storage(),
	storageLastOffset(BITSIZE)
	{}
	void shrink_to_fit()
	{
		storage.shrink_to_fit();
	}
	void reserve(size_t size)
	{
		storage.reserve((size + CHARSINBLOCK - 1) / CHARSINBLOCK);
	}
	size_t size() const
	{
		return storage.size() * CHARSINBLOCK - (BITSIZE - storageLastOffset) / 2;
	}
	void addChar(char c)
	{
		if (storageLastOffset == BITSIZE)
		{
			storage.push_back(0);
			storageLastOffset = 0;
		}
		Word num = 0;
		switch(c)
		{
			case 'a':
			case 'A':
				num = 0;
				break;
			case 'c':
			case 'C':
				num = 1;
				break;
			case 'g':
			case 'G':
				num = 2;
				break;
			case 't':
			case 'T':
				num = 3;
				break;
			default:
				assert(false);
		}
		storage.back() |= num << storageLastOffset;
		storageLastOffset += 2;
	}
	Word getWord(size_t block) const
	{
		return storage[block];
	}
	char getChar(size_t pos) const
	{
		return "ACGT"[getCharI(pos)];
	}
	int getCharI(size_t pos) const
	{
		size_t block = pos / CHARSINBLOCK;
		size_t offset = (pos % CHARSINBLOCK) * 2;
		return (getWord(block) >> offset) & 3;
	}
private:
	std::vector<Word> storage;
	size_t storageLastOffset;
};

#endif
