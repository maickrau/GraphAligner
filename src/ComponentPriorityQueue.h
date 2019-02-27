#ifndef ComponentPriorityQueue_h
#define ComponentPriorityQueue_h

#include <queue>
#include <sparsehash/sparse_hash_map>
#include "ThreadReadAssertion.h"

template <typename T, bool SparseStorage>
class ComponentPriorityQueue
{
	struct PrioritizedItem
	{
		PrioritizedItem(size_t component, int score, size_t index) : component(component), score(score), index(index) {}
		size_t component;
		int score;
		size_t index;
		bool operator<(const PrioritizedItem& other) const { return component < other.component || (component == other.component && score < other.score); }
		bool operator>(const PrioritizedItem& other) const { return component > other.component || (component == other.component && score > other.score); }
	};
public:
	constexpr bool IsComponentPriorityQueue() { return true; }
	ComponentPriorityQueue(size_t maxNode) :
	activeQueues(),
	active(),
	extras()
	{
		initialize(maxNode);
	}
	ComponentPriorityQueue() :
	activeQueues(),
	active(),
	extras()
	{
	}
	template <bool Sparse = SparseStorage>
	typename std::enable_if<Sparse>::type initialize(size_t maxNode)
	{
		extras.set_deleted_key(std::numeric_limits<size_t>::max());
		active.resize(maxNode, false);
	}
	template <bool Sparse = SparseStorage>
	typename std::enable_if<!Sparse>::type initialize(size_t maxNode)
	{
		extras.resize(maxNode);
		active.resize(maxNode, false);
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	T& top()
	{
		assert(activeQueues.size() > 0);
		auto index = activeQueues.top().index;
		assert(active[index]);
		assert(extras[index].size() > 0);
		return extras[index][0];
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	void pop()
	{
		assert(activeQueues.size() > 0);
		size_t index = activeQueues.top().index;
		assert(active[index]);
		assert(extras[index].size() > 0);
		extras[index].clear();
		active[index] = false;
		activeQueues.pop();
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	size_t size() const
	{
		return activeQueues.size();
	}
	void insert(size_t component, const T& item)
	{
		assert(false);
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	void insert(size_t component, int score, const T& item)
	{
		size_t index = getId(item);
		assert(SparseStorage || index < extras.size());
		if (!active[index])
		{
			assert(extras[index].size() == 0);
			activeQueues.emplace(component, score, index);
			active[index] = true;
		}
		extras[index].push_back(item);
	}
	void clear()
	{
		while (activeQueues.size() > 0)
		{
			size_t index = activeQueues.top().index;
			assert(active[index]);
			removeExtras(index);
			active[index] = false;
			activeQueues.pop();
		}
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
	const std::vector<T>& getExtras(size_t index) const
	{
		assert(SparseStorage ||index < extras.size());
		return getVec(extras, index);
	}
	void removeExtras(size_t index)
	{
		assert(SparseStorage ||index < extras.size());
		extras[index].clear();
	}
	size_t extraSize(size_t index) const
	{
		assert(SparseStorage ||index < extras.size());
		return getVec(extras, index).size();
	}
	bool valid() const
	{
		return active.size() > 0;
	}
private:
	const std::vector<T>& getVec(const std::vector<std::vector<T>>& list, size_t index) const
	{
		return list[index];
	}
	const std::vector<T>& getVec(const google::sparse_hash_map<size_t, std::vector<T>>& list, size_t index) const
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
	std::priority_queue<PrioritizedItem, std::vector<PrioritizedItem>, std::greater<PrioritizedItem>> activeQueues;
	std::vector<bool> active;
	typename std::conditional<SparseStorage, google::sparse_hash_map<size_t, std::vector<T>>, std::vector<std::vector<T>>>::type extras;
};

#endif
