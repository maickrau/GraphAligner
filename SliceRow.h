#ifndef SliceRow_h
#define SliceRow_h

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
	size_t size() const
	{
		return numElems;
	}
	size_t numBlocks() const
	{
		return ends.size();
	}
private:
	std::vector<LengthType> ends;
	std::vector<LengthType> starts;
	size_t numElems;
};

#endif
