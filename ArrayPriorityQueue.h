#ifndef ArrayPriorityQueue_h
#define ArrayPriorityQueue_h

#include <queue>
#include "ThreadReadAssertion.h"

template <typename T>
class ArrayPriorityQueue
{
public:
	ArrayPriorityQueue(size_t maxPriority) :
	activeQueues(),
	queues(),
	numItems(0)
	{
		queues.resize(maxPriority);
	}
	T& top()
	{
		assert(activeQueues.size() > 0);
		size_t queue = activeQueues.top();
		assert(queues[queue].size() > 0);
		return queues[queue].back();
	}
	void pop()
	{
		size_t queue = activeQueues.top();
		assert(queues[queue].size() > 0);
		queues[queue].pop_back();
		if (queues[queue].size() == 0) activeQueues.pop();
		numItems--;
	}
	size_t size() const
	{
		return numItems;
	}
	void insert(size_t priority, const T& item)
	{
		assert(priority < queues.size());
		queues[priority].push_back(item);
		if (queues[priority].size() == 1) activeQueues.emplace(priority);
		numItems++;
	}
	void clear()
	{
		while (activeQueues.size() > 0)
		{
			size_t queue = activeQueues.top();
			queues[queue].clear();
			activeQueues.pop();
		}
		numItems = 0;
	}
private:
	std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> activeQueues;
	std::vector<std::vector<T>> queues;
	size_t numItems;
};

#endif
