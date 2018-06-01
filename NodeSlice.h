#ifndef NodeSlice_h
#define NodeSlice_h

#include <memory>
#include <limits>
#include <unordered_map>
#include <vector>
#include "ThreadReadAssertion.h"
#include "WordSlice.h"

struct NodeSliceMapItem
{
	size_t start;
	size_t end;
	int minScore;
	size_t startOffset;
	size_t endOffset;
};

template <typename LengthType, typename ScoreType, typename Word>
class WordContainer
{
public:
	typedef WordSlice<LengthType, ScoreType, Word> Slice;
	class SmallSlice
	{
	public:
		Word VP;
		Word VN;
		uint16_t plusMinScore;
		bool sliceExists;
		int8_t minPlusEndScore;
	};
	class TinySlice
	{
	public:
		uint16_t plusMinScore;
		char VPVNLastBit;
	};
	class WordContainerIterator : std::iterator<std::random_access_iterator_tag, Slice>
	{
	public:
		WordContainerIterator() :
		container(nullptr),
		pos(0)
		{
		}
		WordContainerIterator(WordContainer* container, size_t pos) :
		container(container),
		pos(pos)
		{
		}
		Slice& operator*()
		{
			return (*container)[pos];
		}
		const Slice operator*() const
		{
			return (*container)[pos];
		}
		Slice& operator[](size_t index)
		{
			return (*container)[pos + index];
		}
		const Slice operator[](size_t index) const
		{
			return (*container)[pos + index];
		}
		WordContainerIterator& operator++()
		{
			pos++;
			return *this;
		}
		WordContainerIterator operator++(int)
		{
			auto result = *this;
			pos++;
			return result;
		}
		WordContainerIterator& operator--()
		{
			pos--;
			return *this;
		}
		WordContainerIterator operator--(int)
		{
			auto result = *this;
			pos--;
			return result;
		}
		bool operator==(const WordContainerIterator& other) const
		{
			return container == other.container && pos == other.pos;
		}
		bool operator!=(const WordContainerIterator& other) const
		{
			return !(*this == other);
		}
		bool operator<(const WordContainerIterator& other) const
		{
			return pos < other.pos;
		}
		WordContainerIterator operator+(size_t off) const
		{
			WordContainerIterator result { *this };
			result.pos = pos + off;
			return result;
		}
		WordContainerIterator operator-(size_t off) const
		{
			WordContainerIterator result { *this };
			result.pos = pos - off;
			return result;
		}
		WordContainerIterator& operator+=(size_t off)
		{
			pos += off;
			return *this;
		}
		WordContainerIterator& operator-=(size_t off)
		{
			pos -= off;
			return *this;
		}
		size_t operator-(const WordContainerIterator& other) const
		{
			return pos - other.pos;
		}
		WordContainer* container;
		size_t pos;
	private:
	};
	class WordContainerConstIterator : std::iterator<std::random_access_iterator_tag, Slice>
	{
	public:
		WordContainerConstIterator() :
		container(nullptr),
		pos(0)
		{
		}
		WordContainerConstIterator(const WordContainer* container, size_t pos) :
		container(container),
		pos(pos)
		{
		}
		WordContainerConstIterator(const WordContainerIterator& iter) :
		container(iter.container),
		pos(iter.pos)
		{
		}
		const Slice operator*() const
		{
			return (*container)[pos];
		}
		const Slice operator[](size_t index) const
		{
			return (*container)[pos + index];
		}
		WordContainerConstIterator& operator++()
		{
			pos++;
			return *this;
		}
		WordContainerConstIterator operator++(int)
		{
			auto result = *this;
			pos++;
			return result;
		}
		WordContainerConstIterator& operator--()
		{
			pos--;
			return *this;
		}
		WordContainerConstIterator operator--(int)
		{
			auto result = *this;
			pos--;
			return result;
		}
		bool operator==(const WordContainerConstIterator& other) const
		{
			return container == other.container && pos == other.pos;
		}
		bool operator!=(const WordContainerConstIterator& other) const
		{
			return !(*this == other);
		}
		bool operator<(const WordContainerConstIterator& other) const
		{
			return pos < other.pos;
		}
		WordContainerConstIterator operator+(size_t off) const
		{
			WordContainerConstIterator result { *this };
			result.pos = pos + off;
			return result;
		}
		WordContainerConstIterator operator-(size_t off) const
		{
			WordContainerConstIterator result { *this };
			result.pos = pos - off;
			return result;
		}
		WordContainerConstIterator& operator+=(size_t off)
		{
			pos += off;
			return *this;
		}
		WordContainerConstIterator& operator-=(size_t off)
		{
			pos -= off;
			return *this;
		}
		size_t operator-(const WordContainerConstIterator& other) const
		{
			return pos - other.pos;
		}
	private:
		const WordContainer* container;
		size_t pos;
	};
	using iterator = WordContainerIterator;
	using const_iterator = WordContainerConstIterator;
	class ContainerView
	{
	public:
		ContainerView(iterator begin, iterator end, int minScore) :
		_begin(begin),
		_end(end),
		_cbegin(begin),
		_cend(end),
		_minScore(minScore)
		{
		}
		ContainerView(const_iterator begin, const_iterator end, int minScore) :
		_cbegin(begin),
		_cend(end),
		_minScore(minScore)
		{
		}
		iterator begin()
		{
			return _begin;
		}
		iterator end()
		{
			return _end;
		}
		const_iterator begin() const
		{
			return _cbegin;
		}
		const_iterator end() const
		{
			return _cend;
		}
		Slice& back()
		{
			return *(_end - 1);
		}
		const Slice back() const
		{
			return *(_cend - 1);
		}
		size_t size() const
		{
			return _cend - _cbegin;
		}
		Slice& operator[](size_t pos)
		{
			return *(_begin + pos);
		}
		const Slice operator[](size_t pos) const
		{
			return *(_cbegin + pos);
		}
		int minScore() const
		{
			return _minScore;
		}
	private:
		iterator _begin;
		iterator _end;
		const_iterator _cbegin;
		const_iterator _cend;
		int _minScore;
	};
	WordContainer() :
	mutableSlices()
	{
	}
	ContainerView getView(size_t start, size_t end, int minScore)
	{
		return ContainerView { begin() + start, begin() + end, minScore };
	}
	const ContainerView getView(size_t start, size_t end, int minScore) const
	{
		return ContainerView { begin() + start, begin() + end, minScore };
	}
	Slice& operator[](size_t index)
	{
		return mutableSlices[index];
	}
	const Slice operator[](size_t index) const
	{
		return mutableSlices[index];
	}
	iterator begin()
	{
		return iterator { this, 0 };
	}
	iterator end()
	{
		return iterator { this, size() };
	}
	const_iterator begin() const
	{
		return const_iterator { this, 0 };
	}
	const_iterator end() const
	{
		return const_iterator { this, size() };
	}
	size_t size() const
	{
		return mutableSlices.size();
	}
	void resize(size_t size, Slice value)
	{
		mutableSlices.resize(size, value);
	}
	Slice& back()
	{
		return mutableSlices.back();
	}
	const Slice back() const
	{
		return (*this)[size()-1];
	}
	void reserve(size_t size)
	{
		mutableSlices.reserve(size);
	}
private:
	std::vector<Slice> mutableSlices;
};

template <typename T>
class NodeSlice
{
public:
	using MapType = std::unordered_map<size_t, NodeSliceMapItem>;
	using MapItem = NodeSliceMapItem;
	using Container = WordContainer<size_t, int, uint64_t>;
	using View = Container::ContainerView;
	class NodeSliceIterator : std::iterator<std::forward_iterator_tag, std::pair<size_t, View>>
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
		std::pair<size_t, View> operator*()
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, slice->slices->getView(info.start, info.end, info.minScore));
			}
			else
			{
				return std::make_pair(mappos->first, slice->slices->getView(mappos->second.start, mappos->second.end, mappos->second.minScore));
			}
		}
		const std::pair<size_t, const View> operator*() const
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, slice->slices->getView(info.start, info.end, info.minScore));
			}
			else
			{
				return std::make_pair(mappos->first, slice->slices->getView(mappos->second.start, mappos->second.end, mappos->second.minScore));
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
	class NodeSliceConstIterator : std::iterator<std::forward_iterator_tag, const std::pair<size_t, View>>
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
		const std::pair<size_t, const View> operator*() const
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, slice->slices->getView(info.start, info.end, info.minScore));
			}
			else
			{
				return std::make_pair(mappos->first, slice->slices->getView(mappos->second.start, mappos->second.end, mappos->second.minScore));
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
	nodes(new MapType),
	slices(new Container)
	{
	}
	NodeSlice(std::vector<NodeSliceMapItem>* vectorMap) :
	vectorMap(vectorMap),
	nodes(new MapType),
	slices(new Container)
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
	void reserve(size_t size)
	{
		slices->reserve(size);
	}
	void addNode(size_t nodeIndex, size_t size)
	{
		addNode(nodeIndex, size, typename Container::Slice {0, 0, 0, 0, 0 });
	}
	void addNode(size_t nodeIndex, size_t size, typename Container::Slice value)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert((*vectorMap)[nodeIndex].start == (*vectorMap)[nodeIndex].end);
			(*vectorMap)[nodeIndex] = { slices->size(), slices->size() + size, 0, size, 0 };
			activeVectorMapIndices.push_back(nodeIndex);
		}
		else
		{
			assert(nodes->find(nodeIndex) == nodes->end());
			(*nodes)[nodeIndex] = { slices->size(), slices->size() + size, 0, size, 0 };
		}
		slices->resize(slices->size() + size, value);
	}
	View node(size_t nodeIndex)
	{
		auto item = getMapItem(nodeIndex);
		return slices->getView(item.start, item.end, item.minScore);
	}
	const View node(size_t nodeIndex) const
	{
		auto item = getMapItem(nodeIndex);
		return slices->getView(item.start, item.end, item.minScore);
	}
	bool hasNode(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			return (*vectorMap)[nodeIndex].start != (*vectorMap)[nodeIndex].end;
		}
		else
		{
			return nodes->find(nodeIndex) != nodes->end();
		}
	}
	int minScore(size_t nodeIndex) const
	{
		return getMapItem(nodeIndex).minScore;
	}
	void setMinScoreIfSmaller(size_t nodeIndex, int score)
	{
		auto oldScore = minScore(nodeIndex);
		if (score < oldScore) setMinScore(nodeIndex, score);
	}
	void setMinScore(size_t nodeIndex, int score)
	{
		getMapItem(nodeIndex).minScore = score;
	}
	size_t startIndex(size_t nodeIndex) const
	{
		return getMapItem(nodeIndex).startOffset;
	}
	void setStartIndex(size_t nodeIndex, size_t index)
	{
		getMapItem(nodeIndex).startOffset = index;
	}
	size_t endIndex(size_t nodeIndex) const
	{
		return getMapItem(nodeIndex).endOffset;
	}
	void setEndIndex(size_t nodeIndex, size_t index)
	{
		getMapItem(nodeIndex).endOffset = index;
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
		NodeSliceMapItem empty;
		empty.start = 0;
		empty.end = 0;
		empty.minScore = 0;
		empty.startOffset = 0;
		empty.endOffset = 0;
		for (auto index : activeVectorMapIndices)
		{
			(*vectorMap)[index] = empty;
		}
		activeVectorMapIndices.clear();
	}
	NodeSliceMapItem& getMapItem(size_t nodeIndex)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert((*vectorMap)[nodeIndex].start != (*vectorMap)[nodeIndex].end);
			return (*vectorMap)[nodeIndex];
		}
		else
		{
			auto found = nodes->find(nodeIndex);
			assert(found != nodes->end());
			return found->second;
		}
	}
	const NodeSliceMapItem& getMapItem(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert((*vectorMap)[nodeIndex].start != (*vectorMap)[nodeIndex].end);
			return (*vectorMap)[nodeIndex];
		}
		else
		{
			auto found = nodes->find(nodeIndex);
			assert(found != nodes->end());
			return found->second;
		}
	}
	std::vector<NodeSliceMapItem>* vectorMap;
	std::vector<size_t> activeVectorMapIndices;
	std::shared_ptr<MapType> nodes;
	std::shared_ptr<Container> slices;
	friend class NodeSliceIterator;
	friend class NodeSliceConstIterator;
};

#endif