#ifndef ArrayPriorityQueue_h
#define ArrayPriorityQueue_h

#include <queue>
#include "ThreadReadAssertion.h"

template <typename T>
class ArrayPriorityQueue
{
public:
	ArrayPriorityQueue(size_t maxPriority, size_t maxExtras) :
	activeQueues(),
	extras(),
	queues(),
	numItems(0)
	{
		extras.resize(maxExtras, std::vector<T>{});
		queues.resize(maxPriority);
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	T& top()
	{
		assert(activeQueues.size() > 0);
		size_t queue = activeQueues.top();
		assert(queues[queue].size() > 0);
		return queues[queue].back();
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	void pop()
	{
		size_t queue = activeQueues.top();
		assert(queues[queue].size() > 0);
		queues[queue].pop_back();
		if (queues[queue].size() == 0) activeQueues.pop();
		numItems--;
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	size_t size() const
	{
		return numItems;
	}
	void insert(size_t component, int score, const T& item)
	{
		assert(false);
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	void insert(size_t priority, const T& item)
	{
		assert(priority < queues.size());
		queues[priority].push_back(item);
		assert(getId(item) < extras.size());
		extras[getId(item)].push_back(item);
		if (queues[priority].size() == 1) activeQueues.emplace(priority);
		numItems++;
	}
	void clear()
	{
		while (activeQueues.size() > 0)
		{
			size_t queue = activeQueues.top();
			for (auto item : queues[queue])
			{
				removeExtras(getId(item));
			}
			queues[queue].clear();
			activeQueues.pop();
		}
		numItems = 0;
	}
	const std::vector<T>& getExtras(size_t index) const
	{
		assert(index < extras.size());
		return extras[index];
	}
	void removeExtras(size_t index)
	{
		assert(index < extras.size());
		extras[index].clear();
	}
	size_t extraSize(size_t index) const
	{
		assert(index < extras.size());
		return extras[index].size();
	}
private:
	size_t getId(const T& item) const
	{
		return item.target;
	}
	std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> activeQueues;
	std::vector<std::vector<T>> extras;
	std::vector<std::vector<T>> queues;
	size_t numItems;
};

#endif
