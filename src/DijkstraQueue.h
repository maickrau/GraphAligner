#ifndef DijkstraQueue_h
#define DijkstraQueue_h

#include <queue>
#include <phmap.h>
#include "ThreadReadAssertion.h"

namespace std
{
    template<>
    struct hash<std::pair<size_t, size_t>>
    {
    	size_t operator()(const std::pair<size_t, size_t> pair) const
    	{
    		return std::hash<size_t>{}(pair.first) ^ std::hash<size_t>{}(pair.second);
    	}
    };
}

template <typename T>
class DijkstraPriorityQueue
{
public:
	constexpr bool IsComponentPriorityQueue() { return false; }
	DijkstraPriorityQueue() :
	activeQueues(),
	extras(),
	queues(),
	numItems(0),
	zeroScore(0)
	{
		initialize(129);
	}
	void initialize(size_t maxPriority)
	{
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
		assert(numItems > 0);
		assert(activeQueues.size() > 0);
		size_t queue = activeQueues.top();
		assert(queue < queues.size());
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
		assert(priority >= zeroScore);
		priority -= zeroScore;
		assert(priority < queues.size());
		queues[priority].push_back(item);
		extras[getId(item)].push_back(item);
		if (queues[priority].size() == 1) activeQueues.emplace(priority);
		numItems++;
	}
	void clear()
	{
		for (size_t i = 0; i < queues.size(); i++)
		{
			queues[i].clear();
		}
		typename std::remove_reference<decltype(activeQueues)>::type tmp;
		std::swap(tmp, activeQueues);
		typename std::remove_reference<decltype(extras)>::type tmp2;
		std::swap(tmp2, extras);
		numItems = 0;
		zeroScore = 0;
		sparsify();
	}
	void increaseScore(size_t increase)
	{
		assert(increase > 0);
		assert(increase < queues.size());
		assert(activeQueues.size() == 0 || activeQueues.top() >= increase);
		typename std::remove_reference<decltype(activeQueues)>::type tmp;
		std::swap(tmp, activeQueues);
		for (size_t i = 0; i < queues.size() - increase; i++)
		{
			assert(queues[i].size() == 0);
			std::swap(queues[i], queues[i + increase]);
			if (queues[i].size() > 0) activeQueues.emplace(i);
		}
		for (size_t i = queues.size() - increase; i < queues.size(); i++)
		{
			assert(queues[i].size() == 0);
		}
		zeroScore += increase;
	}
	void sparsify()
	{
		decltype(extras) empty;
		std::swap(extras, empty);
	}
	const std::vector<T>& getExtras(size_t slice, size_t index)
	{
		return getExtras(std::make_pair(slice, index));
	}
	const std::vector<T>& getExtras(std::pair<size_t, size_t> index)
	{
		return getVec(index);
	}
	void removeExtras(size_t slice, size_t index)
	{
		removeExtras(std::make_pair(slice, index));
	}
	void removeExtras(std::pair<size_t, size_t> index)
	{
		extras[index].clear();
	}
	size_t extraSize(size_t slice, size_t index) const
	{
		return extraSize(std::make_pair(slice, index));
	}
	size_t extraSize(std::pair<size_t, size_t> index) const
	{
		return getVec(index).size();
	}
	size_t zero() const
	{
		return zeroScore;
	}
private:
	const std::vector<T>& getVec(std::pair<size_t, size_t> index) const
	{
		static std::vector<T> empty;
		auto found = extras.find(index);
		if (found == extras.end()) return empty;
		return found->second;
	}
	std::pair<size_t, size_t> getId(const T& item) const
	{
		return std::make_pair(item.slice, item.target);
	}
	std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> activeQueues;
	phmap::flat_hash_map<std::pair<size_t, size_t>, std::vector<T>> extras;
	std::vector<std::vector<T>> queues;
	size_t numItems;
	size_t zeroScore;
};

#endif
