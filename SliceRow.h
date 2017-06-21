#ifndef SliceRow_h
#define SliceRow_h

#include "ThreadReadAssertion.h"

template <typename LengthType>
class SliceRow
{
public:
	class SliceRowConstIterator
	{
	public:
		SliceRowConstIterator(const SliceRow& r) : row(r)
		{
		}
		SliceRowConstIterator& operator++()
		{
			if (loc == nextBreak)
			{
				index++;
				assert(index < row.ends.size());
				nextBreak = row.ends[index];
				loc = row.starts[index];
			}
			else
			{
				assert(loc < nextBreak);
				loc++;
			}
			return *this;
		}
		SliceRowConstIterator operator++(int)
		{
			SliceRowConstIterator result = *this;
			(*this)++;
			return result;
		}
		bool operator==(const SliceRowConstIterator& other) const
		{
			return loc == other.loc;
		}
		bool operator!=(const SliceRowConstIterator& other) const
		{
			return loc != other.loc;
		}
		LengthType operator*() const
		{
			return loc;
		}
	private:
		const SliceRow& row;
		LengthType nextBreak;
		LengthType loc;
		size_t index;
		friend class SliceRow;
	};
	typedef SliceRowConstIterator const_iterator;
	SliceRow() : numElems(0) {};
	int count(LengthType val) const
	{
		auto found = std::lower_bound(ends.begin(), ends.end(), val);
		if (found == ends.end()) return 0;
		auto index = found - ends.begin();
		if (starts[index] <= val) return 1;
		return 0;
	}
	void insert(LengthType val)
	{
		//most inserts happen as a contiguous, ordered block, so check this case first
		if (ends.size() > 0 && ends.back() == val-1)
		{
			ends.back()++;
			numElems++;
			return;
		}
		auto found = std::lower_bound(ends.begin(), ends.end(), val);
		//nothing bigger than this exists, but this is next to the last block
		if (found == ends.end() && ends.size() > 0 && ends.back() == val-1)
		{
			ends.back()++;
			numElems++;
			return;
		}
		//nothing bigger than this exists and it's not next to a block
		if (found == ends.end())
		{
			ends.emplace_back(val);
			starts.emplace_back(val);
			numElems++;
			return;
		}
		auto index = found - ends.begin();
		//already exists
		if (starts[index] <= val) return;
		//right after an existing block
		if (index > 0 && ends[index-1] == val-1)
		{
			ends[index-1]++;
			if (ends[index-1] + 1 == starts[index])
			{
				//merge
				starts[index] = starts[index-1];
				ends.erase(ends.begin() + (index-1));
				starts.erase(starts.begin() + (index-1));
			}
			numElems++;
			return;
		}
		//right before an existing block
		if (starts[index] == val + 1)
		{
			starts[index]--;
			if (index > 0 && starts[index] == ends[index-1] + 1)
			{
				//merge
				starts[index] = starts[index-1];
				ends.erase(ends.begin() + (index-1));
				starts.erase(starts.begin() + (index-1));
			}
			numElems++;
			return;
		}
		//between two blocks but not next to either
		ends.insert(ends.begin() + index, val);
		starts.insert(starts.begin() + index, val);
		numElems++;
	}
	const_iterator begin() const
	{
		const_iterator result {*this};
		if (ends.size() == 0)
		{
			result.index = 0;
			result.nextBreak = 0;
			result.loc = 0;
			return result;
		}
		result.index = 0;
		result.nextBreak = ends[0];
		result.loc = starts[0];
		return result;
	}
	const_iterator end() const
	{
		const_iterator result {*this};
		if (ends.size() == 0)
		{
			result.index = 0;
			result.nextBreak = 0;
			result.loc = 0;
			return result;
		}
		result.index = ends.size()-1;
		result.nextBreak = ends[result.index];
		result.loc = ends[result.index];
		return result;
	}
	template <typename Iterator>
	void insert(Iterator begin, Iterator end)
	{
		while (begin != end)
		{
			insert(*begin);
			++begin;
		}
	}
	void insertBlock(LengthType blockStart, LengthType blockEnd)
	{
		assert(blockEnd >= blockStart);
		assert(ends.size() == starts.size());
		//right after the final block
		if (ends.size() > 0 && ends.back() == blockStart-1)
		{
			numElems += blockEnd - ends.back();
			ends.back() = blockEnd;
			return;
		}
		//right before the first block
		if (ends.size() > 0 && starts[0] == blockEnd+1)
		{
			numElems += starts[0] - blockStart;
			starts[0] = blockStart;
			return;
		}
		//after all blocks, not next to the final block
		if (ends.size() > 0 && ends.back() < blockStart-1)
		{
			starts.push_back(blockStart);
			ends.push_back(blockEnd);
			numElems += blockEnd - blockStart + 1;
			return;
		}
		//before all blocks, not next to the first block
		if (ends.size() > 0 && starts[0] > blockEnd+1)
		{
			starts.insert(starts.begin(), blockStart);
			ends.insert(ends.begin(), blockEnd);
			numElems += blockEnd - blockStart + 1;
			return;
		}
		//no existing blocks
		if (ends.size() == 0)
		{
			starts.push_back(blockStart);
			ends.push_back(blockEnd);
			numElems += blockEnd - blockStart + 1;
			return;
		}
		auto foundEnd = std::lower_bound(ends.begin(), ends.end(), blockEnd);
		auto foundStart = std::lower_bound(ends.begin(), ends.end(), blockStart);
		size_t endIndex = foundEnd - ends.begin();
		size_t startIndex = foundStart - ends.begin();
		//todo handle the case where the inserted block overlaps multiple existing blocks
		//overlaps one existing block
		if (endIndex == startIndex+1)
		{
			//only overlaps the block, and doesn't extend before it
			assert(starts[startIndex] <= blockStart);
			//only overlaps exactly one block
			assert(endIndex == ends.size() || starts[endIndex] > blockEnd);
			assert(blockEnd > ends[startIndex]);
			numElems += blockEnd - ends[startIndex];
			ends[startIndex] = blockEnd;
			//merge with next block
			if (endIndex < ends.size() && ends[startIndex]+1 == starts[endIndex])
			{
				starts[endIndex] = starts[startIndex];
				starts.erase(starts.begin() + startIndex);
				ends.erase(ends.begin() + startIndex);
			}
			return;
		}
		assert(endIndex == startIndex);
		//already exists, don't need to do anything
		if (startIndex < ends.size() && starts[startIndex] <= blockStart && ends[startIndex] >= blockEnd)
		{
			return;
		}
		//right before a block
		if (startIndex < ends.size() && ends[startIndex] >= blockEnd && starts[startIndex] <= blockEnd+1)
		{
			assert(startIndex == 0 || ends[startIndex-1] < blockStart);
			assert(blockStart <= starts[startIndex]);
			numElems += starts[startIndex] - blockStart;
			starts[startIndex] = blockStart;
			//merge with the previous block
			if (startIndex > 0 && starts[startIndex] == ends[startIndex-1] + 1)
			{
				ends[startIndex-1] = ends[startIndex];
				ends.erase(ends.begin() + startIndex);
				starts.erase(starts.begin() + startIndex);
			}
			return;
		}
		//right after a block
		if (startIndex > 0 && starts[startIndex-1] <= blockStart && ends[startIndex-1] >= blockStart-1)
		{
			assert(startIndex == endIndex);
			assert(endIndex == ends.size() || starts[endIndex] > blockEnd);
			assert(blockEnd >= ends[startIndex-1]);
			numElems += blockEnd - ends[startIndex-1];
			ends[startIndex-1] = blockEnd;
			//merge with the next block
			if (starts[startIndex] == ends[startIndex-1]+1)
			{
				ends[startIndex-1] = ends[startIndex];
				ends.erase(ends.begin() + startIndex);
				starts.erase(starts.begin() + startIndex);
			}
			return;
		}
		//between two existing blocks, not touching either
		if (startIndex > 0 && startIndex < ends.size() && ends[startIndex-1] < blockStart-1 && starts[startIndex] > blockEnd + 1)
		{
			assert(startIndex == 0 || ends[startIndex-1] < blockStart-1);
			assert(startIndex == ends.size()-1 || starts[startIndex+1] > blockEnd+1);
			starts.insert(starts.begin() + startIndex, blockStart);
			ends.insert(ends.begin() + startIndex, blockEnd);
			numElems += blockEnd - blockStart + 1;
			return;
		}
		//TODO handle the rest of the cases. for now only these are needed.
		assert(false);
	}
	size_t getSolidIndex(size_t column) const
	{
		if (itemsBefore.size() == 0) calcItemsBefore();
		assert(itemsBefore.size() == ends.size());
		auto endIter = std::lower_bound(ends.begin(), ends.end(), column);
		assert(endIter != ends.end());
		size_t index = endIter - ends.begin();
		assert(starts[index] <= column);
		return itemsBefore[index] + column - starts[index];
	}
	size_t size() const
	{
		return numElems;
	}
	size_t numBlocks() const
	{
		return ends.size();
	}
private:
	void calcItemsBefore() const
	{
		assert(starts.size() == ends.size());
		LengthType total = 0;
		itemsBefore.reserve(starts.size());
		for (size_t i = 0; i < starts.size(); i++)
		{
			itemsBefore.push_back(total);
			total += ends[i]-starts[i]+1;
		}
	}
	std::vector<LengthType> ends;
	std::vector<LengthType> starts;
	mutable std::vector<LengthType> itemsBefore;
	size_t numElems;
};

#endif
