#ifndef ComponentPriorityQueue_h
#define ComponentPriorityQueue_h

#include <queue>
#include "ThreadReadAssertion.h"

template <typename T>
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
	ComponentPriorityQueue(size_t maxNode) :
	activeQueues(),
	active(),
	extras()
	{
		active.resize(maxNode, false);
		extras.resize(maxNode);
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
		assert(index < extras.size());
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
	bool valid() const
	{
		return active.size() > 0;
	}
private:
	size_t getId(const T& item) const
	{
		return item.target;
	}
	std::priority_queue<PrioritizedItem, std::vector<PrioritizedItem>, std::greater<PrioritizedItem>> activeQueues;
	std::vector<bool> active;
	std::vector<std::vector<T>> extras;
};

#endif
