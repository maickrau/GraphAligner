#ifndef NodeSlice_h
#define NodeSlice_h

#include <memory>
#include <limits>
#include <unordered_map>
#include <vector>
#include <sparsehash/dense_hash_map>
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
	using MapType = google::dense_hash_map<size_t, NodeSliceMapItem>;
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
	nodes(nullptr)
	{
	}
	NodeSlice(std::vector<NodeSliceMapItem>* vectorMap) :
	vectorMap(vectorMap),
	nodes(nullptr)
	{
	}
	void addEmptyNodeMap(size_t size)
	{
		assert(nodes == nullptr);
		nodes = std::make_shared<MapType>();
		nodes->set_empty_key(-1);
		nodes->resize(size);
	}
	NodeSlice getMapSlice() const
	{
		if (vectorMap == nullptr) return *this;
		assert(nodes != nullptr);
		assert(nodes->size() == activeVectorMapIndices.size());
		NodeSlice result;
		result.nodes = nodes;
		return result;
	}
	void createNodeMap()
	{
		assert(vectorMap != nullptr);
		addEmptyNodeMap(activeVectorMapIndices.size());
		for (auto index : activeVectorMapIndices)
		{
			assert((*vectorMap)[index].exists);
			(*nodes)[index] = (*vectorMap)[index];
		}
	}
	void removeVectorArray()
	{
		if (vectorMap == nullptr) return;
		assert(nodes != nullptr);
		assert(vectorMap != nullptr);
		assert(nodes->size() == activeVectorMapIndices.size());
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
		assert(nodes != nullptr);
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
			assert(nodes != nullptr);
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
			assert(nodes != nullptr);
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
			assert(nodes != nullptr);
			auto found = nodes->find(nodeIndex);
			if (found == nodes->end()) return false;
			assert(found->second.exists);
			return true;
		}
	}
	void removeNonExistant()
	{
		assert(vectorMap != nullptr);
		std::vector<size_t> newActiveVectorMapIndices;
		newActiveVectorMapIndices.reserve(activeVectorMapIndices.size());
		for (auto index : activeVectorMapIndices)
		{
			if ((*vectorMap)[index].exists) newActiveVectorMapIndices.push_back(index);
		}
		activeVectorMapIndices = newActiveVectorMapIndices;
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
			return NodeSliceIterator { this, emptyIterator, 0 };
		}
		else
		{
			assert(nodes != nullptr);
			return NodeSliceIterator { this, nodes->begin() };
		}
	}
	NodeSliceIterator end()
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceIterator { this, emptyIterator, activeVectorMapIndices.size() };
		}
		else
		{
			assert(nodes != nullptr);
			return NodeSliceIterator { this, nodes->end() };
		}
	}
	NodeSliceConstIterator begin() const
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceConstIterator { this, emptyIterator, 0 };
		}
		else
		{
			assert(nodes != nullptr);
			return NodeSliceConstIterator { this, nodes->begin() };
		}
	}
	NodeSliceConstIterator end() const
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceConstIterator { this, emptyIterator, activeVectorMapIndices.size() };
		}
		else
		{
			assert(nodes != nullptr);
			return NodeSliceConstIterator { this, nodes->end() };
		}
	}
private:
	static typename MapType::iterator emptyIterator;
	void clearVectorMap()
	{
		assert(nodes != nullptr);
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

template <typename LengthType, typename ScoreType, typename Word>
typename NodeSlice<LengthType, ScoreType, Word>::MapType::iterator NodeSlice<LengthType, ScoreType, Word>::emptyIterator = NodeSlice<LengthType, ScoreType, Word>::MapType{}.end();

#endif