#ifndef NodeSlice_h
#define NodeSlice_h

#include <memory>
#include <limits>
#include <unordered_map>
#include <vector>
#include "AlignmentGraph.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"

template <typename LengthType, typename ScoreType, typename Word>
class NodeSlice
{
public:
	struct NodeSliceMapItem
	{
		static constexpr size_t NUM_CHUNKS = (AlignmentGraph::SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
		NodeSliceMapItem() :
		startSlice(),
		endSlice(),
		exists(false),
		HP(),
		HN(),
		minScore(0)
		{
		}
		WordSlice<LengthType, ScoreType, Word> startSlice;
		WordSlice<LengthType, ScoreType, Word> endSlice;
		bool exists;
		Word HP[NUM_CHUNKS];
		Word HN[NUM_CHUNKS];
		ScoreType minScore;
	};
	using MapType = std::unordered_map<size_t, NodeSliceMapItem>;
	using MapItem = NodeSliceMapItem;
	class NodeSliceIterator : std::iterator<std::forward_iterator_tag, std::pair<size_t, MapItem>>
	{
		using map_iterator = typename MapType::iterator;
	public:
		NodeSliceIterator(NodeSlice* slice, map_iterator pos) :
		slice(slice),
		mappos(pos),
		indexPos(-1)
		{
		}
		NodeSliceIterator(NodeSlice* slice, map_iterator pos, size_t indexPos) :
		slice(slice),
		mappos(pos),
		indexPos(indexPos)
		{
		}
		std::pair<size_t, MapItem> operator*()
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, info);
			}
			else
			{
				return std::make_pair(mappos->first, mappos->second);
			}
		}
		const std::pair<size_t, const MapItem> operator*() const
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, info);
			}
			else
			{
				return std::make_pair(mappos->first, mappos->second);
			}
		}
		NodeSliceIterator& operator++()
		{
			if (indexPos != -1)
			{
				indexPos++;
			}
			else
			{
				++mappos;
			}
			return *this;
		}
		bool operator==(const NodeSliceIterator& other) const
		{
			return slice == other.slice && mappos == other.mappos && indexPos == other.indexPos;
		}
		bool operator!=(const NodeSliceIterator& other) const
		{
			return !(*this == other);
		}
	private:
		NodeSlice* slice;
		map_iterator mappos;
		size_t indexPos;
	};
	class NodeSliceConstIterator : std::iterator<std::forward_iterator_tag, const std::pair<size_t, MapItem>>
	{
		using map_iterator = typename MapType::const_iterator;
	public:
		NodeSliceConstIterator(const NodeSlice* slice, map_iterator pos) :
		slice(slice),
		mappos(pos),
		indexPos(-1)
		{
		}
		NodeSliceConstIterator(const NodeSlice* slice, map_iterator pos, size_t indexPos) :
		slice(slice),
		mappos(pos),
		indexPos(indexPos)
		{
		}
		const std::pair<size_t, const MapItem> operator*() const
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, info);
			}
			else
			{
				return std::make_pair(mappos->first, mappos->second);
			}
		}
		NodeSliceConstIterator& operator++()
		{
			if (indexPos != -1)
			{
				indexPos++;
			}
			else
			{
				++mappos;
			}
			return *this;
		}
		bool operator==(const NodeSliceConstIterator& other) const
		{
			return slice == other.slice && mappos == other.mappos && indexPos == other.indexPos;
		}
		bool operator!=(const NodeSliceConstIterator& other) const
		{
			return !(*this == other);
		}
	private:
		const NodeSlice* slice;
		map_iterator mappos;
		size_t indexPos;
	};
	NodeSlice() :
	vectorMap(nullptr),
	nodes(new MapType)
	{
	}
	NodeSlice(std::vector<NodeSliceMapItem>* vectorMap) :
	vectorMap(vectorMap),
	nodes(new MapType)
	{
	}
	void convertVectorArrayToMap()
	{
		assert(nodes->size() == 0);
		assert(vectorMap != nullptr);
		nodes->reserve(activeVectorMapIndices.size());
		for (auto index : activeVectorMapIndices)
		{
			(*nodes)[index] = (*vectorMap)[index];
		}
		clearVectorMap();
		vectorMap = nullptr;
	}
	void addNode(size_t nodeIndex)
	{
		assert(vectorMap != nullptr);
		assert(nodeIndex < vectorMap->size());
		assert(!(*vectorMap)[nodeIndex].exists);
		activeVectorMapIndices.push_back(nodeIndex);
		(*vectorMap)[nodeIndex].minScore = std::numeric_limits<ScoreType>::max();
		(*vectorMap)[nodeIndex].startSlice = { 0, 0, std::numeric_limits<ScoreType>::max() };
		(*vectorMap)[nodeIndex].endSlice = { 0, 0, std::numeric_limits<ScoreType>::max() };
	}
	void addNodeToMap(size_t nodeIndex)
	{
		assert(vectorMap == nullptr);
		(*nodes)[nodeIndex] = {};
	}
	NodeSliceMapItem& node(size_t nodeIndex)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			return (*vectorMap)[nodeIndex];
		}
		else
		{
			auto found = nodes->find(nodeIndex);
			assert(found != nodes->end());
			return found->second;
		}
	}
	const NodeSliceMapItem& node(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			return (*vectorMap)[nodeIndex];
		}
		else
		{
			auto found = nodes->find(nodeIndex);
			assert(found != nodes->end());
			return found->second;
		}
	}
	bool hasNode(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			return (*vectorMap)[nodeIndex].exists;
		}
		else
		{
			auto found = nodes->find(nodeIndex);
			if (found == nodes->end()) return false;
			assert(found->second.exists);
			return true;
		}
	}
	int minScore(size_t nodeIndex) const
	{
		return node(nodeIndex).minScore;
	}
	void setMinScoreIfSmaller(size_t nodeIndex, int score)
	{
		auto oldScore = minScore(nodeIndex);
		if (score < oldScore) setMinScore(nodeIndex, score);
	}
	void setMinScore(size_t nodeIndex, int score)
	{
		node(nodeIndex).minScore = score;
	}
	size_t size() const
	{
		if (vectorMap != nullptr)
		{
			return activeVectorMapIndices.size();
		}
		else
		{
			return nodes->size();
		}
	}
	NodeSliceIterator begin()
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceIterator { this, nodes->end(), 0 };
		}
		else
		{
			return NodeSliceIterator { this, nodes->begin() };
		}
	}
	NodeSliceIterator end()
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceIterator { this, nodes->end(), activeVectorMapIndices.size() };
		}
		else
		{
			return NodeSliceIterator { this, nodes->end() };
		}
	}
	NodeSliceConstIterator begin() const
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceConstIterator { this, nodes->end(), 0 };
		}
		else
		{
			return NodeSliceConstIterator { this, nodes->begin() };
		}
	}
	NodeSliceConstIterator end() const
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceConstIterator { this, nodes->end(), activeVectorMapIndices.size() };
		}
		else
		{
			return NodeSliceConstIterator { this, nodes->end() };
		}
	}
private:
	void clearVectorMap()
	{
		assert(vectorMap != nullptr);
		for (auto index : activeVectorMapIndices)
		{
			(*vectorMap)[index].exists = false;
		}
		activeVectorMapIndices.clear();
	}
	std::vector<NodeSliceMapItem>* vectorMap;
	std::vector<size_t> activeVectorMapIndices;
	std::shared_ptr<MapType> nodes;
	friend class NodeSliceIterator;
	friend class NodeSliceConstIterator;
};

#endif