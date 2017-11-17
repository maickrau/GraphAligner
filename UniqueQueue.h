#ifndef UniqueQueue_h
#define UniqueQueue_h

#include <vector>

template <typename T>
class UniqueQueue
{
public:
	UniqueQueue(size_t maxsize)
	{
		hasItem.resize(maxsize, false);
	}
	void insert(T item)
	{
		if (hasItem[item]) return;
		hasItem[item] = true;
		items.push_back(item);
	}
	size_t size() const
	{
		return items.size();
	}
	void pop()
	{
		assert(hasItem[items.back()]);
		hasItem[items.back()] = false;
		items.pop_back();
	}
	T top() const
	{
		return items.back();
	}
	T operator[](size_t pos) const
	{
		return items[pos];
	}
	void clear()
	{
		while (size() > 0) pop();
	}
	auto begin() const
	{
		return items.begin();
	}
	auto end() const
	{
		return items.end();
	}
	auto begin()
	{
		return items.begin();
	}
	auto end()
	{
		return items.end();
	}
	template <typename Iterator>
	void insert(Iterator start, Iterator end)
	{
		for (; start != end; ++start)
		{
			insert(*start);
		}
	}
private:
	std::vector<T> items;
	std::vector<bool> hasItem;
};

#endif
