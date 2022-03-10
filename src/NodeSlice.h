#ifndef NodeSlice_h
#define NodeSlice_h

#include <memory>
#include <limits>
#include <vector>
#include <type_traits>
#include <phmap.h>
#include "AlignmentGraph.h"
#include "ThreadReadAssertion.h"
#include "WordSlice.h"


template <typename LengthType, typename ScoreType, typename Word>
struct NodeSliceMapItemStruct
{
	static constexpr size_t NUM_CHUNKS = (AlignmentGraph::SPLIT_NODE_SIZE + WordConfiguration<Word>::WordSize - 1) / WordConfiguration<Word>::WordSize;
	NodeSliceMapItemStruct() :
	startSlice(),
	endSlice(),
	exists(false),
	HP(),
	HN(),
	minScore(0)
#ifdef SLICEVERBOSE
	,firstSlicesCalcedWhenCalced(std::numeric_limits<size_t>::max())
	,slicesCalcedWhenCalced(std::numeric_limits<size_t>::max())
#endif
	{
		for (size_t i = 0; i < NUM_CHUNKS; i++)
		{
			HP[i] = 0;
			HN[i] = 0;
		}
	}
	WordSlice<LengthType, ScoreType, Word> startSlice;
	WordSlice<LengthType, ScoreType, Word> endSlice;
	bool exists;
	Word HP[NUM_CHUNKS];
	Word HN[NUM_CHUNKS];
	ScoreType minScore;
#ifdef SLICEVERBOSE
	size_t firstSlicesCalcedWhenCalced;
	size_t slicesCalcedWhenCalced;
#endif
};

template <typename LengthType, typename ScoreType, typename Word, bool UseVectorMap>
class NodeSlice
{
public:
	using NodeSliceMapItem = NodeSliceMapItemStruct<LengthType, ScoreType, Word>;
	using MapType = phmap::flat_hash_map<size_t, NodeSliceMapItem>;
	using MapItem = NodeSliceMapItem;
	class NodeSliceIterator : std::iterator<std::forward_iterator_tag, std::pair<size_t, MapItem>>
	{
		using map_iterator = typename MapType::iterator;
	public:
		template <bool HasVectorMap = UseVectorMap>
		NodeSliceIterator(NodeSlice* slice, typename std::enable_if<!HasVectorMap, map_iterator>::type pos) :
		slice(slice),
		mappos(pos)
		{
		}
		template <bool HasVectorMap = UseVectorMap>
		NodeSliceIterator(NodeSlice* slice, typename std::enable_if<HasVectorMap, size_t>::type indexPos) :
		slice(slice),
		indexPos(indexPos)
		{
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, std::pair<size_t, MapItem>>::type operator*()
		{
			auto nodeindex = slice->activeVectorMapIndices[indexPos];
			auto info = (*(slice->vectorMap))[nodeindex];
			return std::make_pair(nodeindex, info);
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, std::pair<size_t, MapItem>>::type operator*()
		{
			return std::make_pair(mappos->first, mappos->second);
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, std::pair<size_t, const MapItem>>::type operator*() const
		{
			auto nodeindex = slice->activeVectorMapIndices[indexPos];
			auto info = (*(slice->vectorMap))[nodeindex];
			return std::make_pair(nodeindex, info);
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, std::pair<size_t, const MapItem>>::type operator*() const
		{
			return std::make_pair(mappos->first, mappos->second);
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, NodeSliceIterator&>::type operator++()
		{
			indexPos++;
			return *this;
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, NodeSliceIterator&>::type operator++()
		{
			++mappos;
			return *this;
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, bool>::type operator==(const NodeSliceIterator& other) const
		{
			assert(slice == other.slice);
			return indexPos == other.indexPos;
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, bool>::type operator==(const NodeSliceIterator& other) const
		{
			assert(slice == other.slice);
			return mappos == other.mappos;
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
		template <bool HasVectorMap = UseVectorMap>
		NodeSliceConstIterator(const NodeSlice* slice, typename std::enable_if<!HasVectorMap, map_iterator>::type pos) :
		slice(slice),
		mappos(pos)
		{
		}
		template <bool HasVectorMap = UseVectorMap>
		NodeSliceConstIterator(const NodeSlice* slice, typename std::enable_if<HasVectorMap, size_t>::type indexPos) :
		slice(slice),
		indexPos(indexPos)
		{
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, const std::pair<size_t, const MapItem>>::type operator*() const
		{
			auto nodeindex = slice->activeVectorMapIndices[indexPos];
			auto info = (*(slice->vectorMap))[nodeindex];
			return std::make_pair(nodeindex, info);
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, const std::pair<size_t, const MapItem>>::type operator*() const
		{
			return std::make_pair(mappos->first, mappos->second);
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, NodeSliceConstIterator&>::type operator++()
		{
			indexPos++;
			return *this;
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, NodeSliceConstIterator&>::type operator++()
		{
			++mappos;
			return *this;
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<HasVectorMap, bool>::type operator==(const NodeSliceConstIterator& other) const
		{
			assert(slice == other.slice);
			return indexPos == other.indexPos;
		}
		template <bool HasVectorMap = UseVectorMap>
		typename std::enable_if<!HasVectorMap, bool>::type operator==(const NodeSliceConstIterator& other) const
		{
			assert(slice == other.slice);
			return mappos == other.mappos;
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
	template <bool HasVectorMap = UseVectorMap>
	NodeSlice(typename std::enable_if<HasVectorMap, std::vector<NodeSliceMapItem>*>::type vectorMap) :
	vectorMap(vectorMap),
	nodes(nullptr)
	{
	}
	void addEmptyNodeMap(size_t size)
	{
		assert(nodes == nullptr);
		nodes = std::make_shared<MapType>();
		nodes->reserve(size);
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, NodeSlice<LengthType, ScoreType, Word, false>>::type getMapSlice() const
	{
		assert(vectorMap != nullptr);
		NodeSlice<LengthType, ScoreType, Word, false> result;
		result.addEmptyNodeMap(activeVectorMapIndices.size());
		for (auto index : activeVectorMapIndices)
		{
			assert((*vectorMap)[index].exists);
			(*result.nodes)[index] = (*vectorMap)[index];
		}
		return result;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap>::type removeVectorArray()
	{
		if (vectorMap == nullptr) return;
		assert(vectorMap != nullptr);
		clearVectorMap();
		vectorMap = nullptr;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap>::type addNode(size_t nodeIndex)
	{
		assert(vectorMap != nullptr);
		assert(nodeIndex < vectorMap->size());
		assert(!(*vectorMap)[nodeIndex].exists);
		activeVectorMapIndices.push_back(nodeIndex);
		(*vectorMap)[nodeIndex].minScore = std::numeric_limits<ScoreType>::max();
		(*vectorMap)[nodeIndex].startSlice = { 0, 0, std::numeric_limits<ScoreType>::max() };
		(*vectorMap)[nodeIndex].endSlice = { 0, 0, std::numeric_limits<ScoreType>::max() };
#ifdef SLICEVERBOSE
		(*vectorMap)[nodeIndex].slicesCalcedWhenCalced = std::numeric_limits<size_t>::max();
		(*vectorMap)[nodeIndex].firstSlicesCalcedWhenCalced = std::numeric_limits<size_t>::max();
#endif
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap>::type addNode(size_t nodeIndex)
	{
		addNodeToMap(nodeIndex);
	}
	void addNodeToMap(size_t nodeIndex)
	{
		assert(nodes != nullptr);
		assert(vectorMap == nullptr);
		auto& node = (*nodes)[nodeIndex];
		node.minScore = std::numeric_limits<ScoreType>::max();
		node.startSlice = { 0, 0, std::numeric_limits<ScoreType>::max() };
		node.endSlice = { 0, 0, std::numeric_limits<ScoreType>::max() };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, NodeSliceMapItem&>::type node(size_t nodeIndex)
	{
		assert(vectorMap != nullptr);
		assert(nodeIndex < vectorMap->size());
		return (*vectorMap)[nodeIndex];
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, NodeSliceMapItem&>::type node(size_t nodeIndex)
	{
		assert(nodes != nullptr);
		auto found = nodes->find(nodeIndex);
		assert(found != nodes->end());
		return found->second;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, const NodeSliceMapItem&>::type node(size_t nodeIndex) const
	{
		assert(vectorMap != nullptr);
		assert(nodeIndex < vectorMap->size());
		return (*vectorMap)[nodeIndex];
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, const NodeSliceMapItem&>::type node(size_t nodeIndex) const
	{
		assert(nodes != nullptr);
		auto found = nodes->find(nodeIndex);
		assert(found != nodes->end());
		return found->second;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, bool>::type hasNode(size_t nodeIndex) const
	{
		assert(vectorMap != nullptr);
		assert(nodeIndex < vectorMap->size());
		return (*vectorMap)[nodeIndex].exists;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, bool>::type hasNode(size_t nodeIndex) const
	{
		assert(nodes != nullptr);
		auto found = nodes->find(nodeIndex);
		if (found == nodes->end()) return false;
		assert(found->second.exists);
		return true;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap>::type removeNonExistant()
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
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap>::type removeNonExistant()
	{
		assert(nodes != nullptr);
		std::vector<std::pair<size_t, MapItem>> newActiveMapItems;
		newActiveMapItems.reserve(activeVectorMapIndices.size());
		for (auto item : (*nodes))
		{
			if (item.second.exists) newActiveMapItems.push_back(item);
		}
		nodes = std::make_shared<MapType>();
		nodes->resize(newActiveMapItems.size());
		for (auto item : newActiveMapItems)
		{
			(*nodes)[item.first] = item.second;
		}
	}
	ScoreType minScore(size_t nodeIndex) const
	{
		return node(nodeIndex).minScore;
	}
	void setMinScoreIfSmaller(size_t nodeIndex, ScoreType score)
	{
		auto oldScore = minScore(nodeIndex);
		if (score < oldScore) setMinScore(nodeIndex, score);
	}
	void setMinScore(size_t nodeIndex, ScoreType score)
	{
		node(nodeIndex).minScore = score;
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, size_t>::type size() const
	{
		assert(vectorMap != nullptr);
		return activeVectorMapIndices.size();
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, size_t>::type size() const
	{
		assert(nodes != nullptr);
		return nodes->size();
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, NodeSliceIterator>::type begin()
	{
		assert(vectorMap != nullptr);
		return NodeSliceIterator { this, 0 };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, NodeSliceIterator>::type begin()
	{
		assert(nodes != nullptr);
		return NodeSliceIterator { this, nodes->begin() };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, NodeSliceIterator>::type end()
	{
		assert(vectorMap != nullptr);
		return NodeSliceIterator { this, activeVectorMapIndices.size() };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, NodeSliceIterator>::type end()
	{
		assert(nodes != nullptr);
		return NodeSliceIterator { this, nodes->end() };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, NodeSliceConstIterator>::type begin() const
	{
		assert(vectorMap != nullptr);
		return NodeSliceConstIterator { this, 0 };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, NodeSliceConstIterator>::type begin() const
	{
		assert(nodes != nullptr);
		return NodeSliceConstIterator { this, nodes->begin() };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap, NodeSliceConstIterator>::type end() const
	{
		assert(vectorMap != nullptr);
		return NodeSliceConstIterator { this, activeVectorMapIndices.size() };
	}
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<!HasVectorMap, NodeSliceConstIterator>::type end() const
	{
		assert(nodes != nullptr);
		return NodeSliceConstIterator { this, nodes->end() };
	}
	bool hasVectorMapCurrently() const
	{
		return vectorMap != nullptr;
	}
	void clear()
	{
		assert(!hasVectorMapCurrently());
		assert(activeVectorMapIndices.size() == 0);
		if (nodes != nullptr)
		{
			nodes->clear();
		}
		else
		{
			addEmptyNodeMap(1);
		}
	}
private:
	template <bool HasVectorMap = UseVectorMap>
	typename std::enable_if<HasVectorMap>::type clearVectorMap()
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
	friend class NodeSlice<LengthType, ScoreType, Word, true>;
};

#endif