#ifndef ArrayPriorityQueue_h
#define ArrayPriorityQueue_h

#include <queue>
#include <phmap.h>
#include "ThreadReadAssertion.h"

template <typename T, bool SparseStorage>
class ArrayPriorityQueue
{
public:
	constexpr bool IsComponentPriorityQueue() { return false; }
	ArrayPriorityQueue(size_t maxPriority, size_t maxExtras) :
	activeQueues(),
	extras(),
	queues(),
	numItems(0)
	{
		initialize(maxPriority, maxExtras);
	}
	ArrayPriorityQueue() :
	activeQueues(),
	extras(),
	queues(),
	numItems(0)
	{
	}
	template <bool Sparse = SparseStorage>
	typename std::enable_if<Sparse>::type initialize(size_t maxPriority, size_t maxExtras)
	{
		queues.resize(maxPriority);
	}
	template <bool Sparse = SparseStorage>
	typename std::enable_if<!Sparse>::type initialize(size_t maxPriority, size_t maxExtras)
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
	void insert(size_t component, int64_t score, const T& item)
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
		assert(SparseStorage || getId(item) < extras.size());
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
		sparsify();
	}

	template<bool Sparse = SparseStorage>
	typename std::enable_if<Sparse>::type sparsify()
	{
		decltype(extras) empty;
		std::swap(extras, empty);
	}
	template<bool Sparse = SparseStorage>
	typename std::enable_if<!Sparse>::type sparsify()
	{
	}
	const std::vector<T>& getExtras(size_t index)
	{
		assert(SparseStorage || index < extras.size());
		return getVec(extras, index);
	}
	void removeExtras(size_t index)
	{
		assert(SparseStorage || index < extras.size());
		extras[index].clear();
	}
	size_t extraSize(size_t index) const
	{
		assert(SparseStorage || index < extras.size());
		return getVec(extras, index).size();
	}
private:
	const std::vector<T>& getVec(const std::vector<std::vector<T>>& list, size_t index) const
	{
		return list[index];
	}
	const std::vector<T>& getVec(const phmap::flat_hash_map<size_t, std::vector<T>>& list, size_t index) const
	{
		static std::vector<T> empty;
		auto found = list.find(index);
		if (found == list.end()) return empty;
		return found->second;
	}
	size_t getId(const T& item) const
	{
		return item.target;
	}
	std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> activeQueues;
	typename std::conditional<SparseStorage, phmap::flat_hash_map<size_t, std::vector<T>>, std::vector<std::vector<T>>>::type extras;
	std::vector<std::vector<T>> queues;
	size_t numItems;
};

#endif
