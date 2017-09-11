#ifndef NodeSlice_h
#define NodeSlice_h

#include <unordered_map>
#include <vector>
#include "ThreadReadAssertion.h"

template <typename T>
class NodeSlice
{
public:
	void addNode(size_t nodeIndex, size_t size)
	{
		nodes[nodeIndex].resize(size);
	}
	std::vector<T>& node(size_t nodeIndex)
	{
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return found->second;
	}
	const std::vector<T>& node(size_t nodeIndex) const
	{
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return found->second;
	}
	auto begin()
	{
		return nodes.begin();
	}
	auto end()
	{
		return nodes.end();
	}
private:
	std::unordered_map<size_t, std::vector<T>> nodes;
};

#endif