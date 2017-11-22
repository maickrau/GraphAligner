#ifndef NodeSlice_h
#define NodeSlice_h

#include <limits>
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
	bool hasNode(size_t nodeIndex) const
	{
		return nodes.find(nodeIndex) != nodes.end();
	}
	int minScore(size_t nodeIndex) const
	{
		auto found = minScores.find(nodeIndex);
		if (found == minScores.end()) return std::numeric_limits<int>::max();
		return found->second;
	}
	void setMinScore(size_t nodeIndex, int score)
	{
		minScores[nodeIndex] = score;
	}
	size_t size() const
	{
		return nodes.size();
	}
	auto begin()
	{
		return nodes.begin();
	}
	auto end()
	{
		return nodes.end();
	}
	auto begin() const
	{
		return nodes.begin();
	}
	auto end() const
	{
		return nodes.end();
	}
private:
	std::unordered_map<size_t, std::vector<T>> nodes;
	std::unordered_map<size_t, int> minScores;
};

#endif