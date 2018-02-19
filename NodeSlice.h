#ifndef NodeSlice_h
#define NodeSlice_h

#include <limits>
#include <unordered_map>
#include <vector>
#include "ThreadReadAssertion.h"
#include "WordSlice.h"

template <typename LengthType, typename ScoreType, typename Word>
class WordContainer
{
public:
	enum FrozenType { FullyMutable, FrozenScores, FrozenLastRow };
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
	minEndScore(0),
	minStartScore(0),
	frozen(FullyMutable)
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
		assert(frozen == FullyMutable);
		return mutableSlices[index];
	}
	const Slice operator[](size_t index) const
	{
		return getSlice(index);
	}
	const Slice getSlice(size_t index) const
	{
		if (frozen == FrozenScores)
		{
			Slice result { frozenSlices[index].VP, frozenSlices[index].VN, 0, minStartScore + frozenSlices[index].plusMinScore, false };
			result.scoreEnd = result.scoreBeforeStart + frozenSlices[index].minPlusEndScore;
			result.sliceExists = frozenSlices[index].sliceExists;
			return result;
		}
		else if (frozen == FrozenLastRow)
		{
			bool VP = frozenSqrtSlices[index].VPVNLastBit & 1;
			bool VN = frozenSqrtSlices[index].VPVNLastBit & 2;
			Slice result { (Word)VP << 63, (Word)VN << 63, minEndScore + frozenSqrtSlices[index].plusMinScore, minEndScore + frozenSqrtSlices[index].plusMinScore - (VP ? 1 : 0) + (VN ? 1 : 0), false };
			if (frozenSqrtSlices[index].VPVNLastBit & 4)
			{
				result.sliceExists = true;
			}
			return result;
		}
		else
		{
			return mutableSlices[index];
		}
	}
	WordContainer getFrozenScores() const
	{
		if (frozen == FrozenScores) return *this;
		assert(frozen == FullyMutable);
		WordContainer result;
		result.frozen = FrozenScores;
		result.frozenSlices.resize(mutableSlices.size());
		result.minStartScore = mutableSlices[0].scoreBeforeStart;
		for (size_t i = 1; i < mutableSlices.size(); i++)
		{
			result.minStartScore = std::min(result.minStartScore, mutableSlices[i].scoreBeforeStart);
		}
		for (size_t i = 0; i < mutableSlices.size(); i++)
		{
			result.frozenSlices[i].VP = mutableSlices[i].VP;
			result.frozenSlices[i].VN = mutableSlices[i].VN;
			assert(mutableSlices[i].scoreBeforeStart >= result.minStartScore);
			assert(mutableSlices[i].scoreBeforeStart - result.minStartScore < std::numeric_limits<decltype(frozenSlices[i].plusMinScore)>::max());
			result.frozenSlices[i].plusMinScore = mutableSlices[i].scoreBeforeStart - result.minStartScore;
			result.frozenSlices[i].minPlusEndScore = mutableSlices[i].scoreEnd - mutableSlices[i].scoreBeforeStart;
			result.frozenSlices[i].sliceExists = mutableSlices[i].sliceExists;
		}
		return result;
	}
	WordContainer getFrozenSqrtEndScores() const
	{
		if (frozen == FrozenLastRow) return *this;
		assert(frozen == FullyMutable);
		WordContainer result;
		result.frozen = FrozenLastRow;
		result.frozenSqrtSlices.resize(mutableSlices.size());
		result.minEndScore = mutableSlices[0].scoreEnd;
		for (size_t i = 1; i < mutableSlices.size(); i++)
		{
			result.minEndScore = std::min(result.minEndScore, mutableSlices[i].scoreEnd);
		}
		for (size_t i = 0; i < mutableSlices.size(); i++)
		{
			result.frozenSqrtSlices[i].VPVNLastBit = 0;
			result.frozenSqrtSlices[i].VPVNLastBit |= mutableSlices[i].VP >> 63;
			result.frozenSqrtSlices[i].VPVNLastBit |= (mutableSlices[i].VN >> 62) & 2;
			result.frozenSqrtSlices[i].VPVNLastBit |= (mutableSlices[i].sliceExists) ? 4 : 0;
			assert(mutableSlices[i].scoreEnd >= result.minEndScore);
			assert(mutableSlices[i].scoreEnd - result.minEndScore < std::numeric_limits<decltype(frozenSqrtSlices[i].plusMinScore)>::max());
			result.frozenSqrtSlices[i].plusMinScore = mutableSlices[i].scoreEnd - result.minEndScore;
		}
		return result;
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
		return frozen == FullyMutable ? mutableSlices.size() : frozen == FrozenScores ? frozenSlices.size() : frozenSqrtSlices.size();
	}
	void resize(size_t size)
	{
		assert(frozen == FullyMutable);
		mutableSlices.resize(size);
	}
	Slice& back()
	{
		assert(frozen == FullyMutable);
		return mutableSlices.back();
	}
	const Slice back() const
	{
		return (*this)[size()-1];
	}
	void reserve(size_t size)
	{
		assert(frozen == FullyMutable);
		mutableSlices.reserve(size);
	}
	ScoreType minScore;
private:
	ScoreType minEndScore;
	ScoreType minStartScore;
	FrozenType frozen;
	std::vector<Slice> mutableSlices;
	std::vector<SmallSlice> frozenSlices;
	std::vector<TinySlice> frozenSqrtSlices;
};

template <typename T>
class NodeSlice
{
public:
	using MapItem = std::tuple<size_t, size_t, int, size_t>;
	using Container = WordContainer<size_t, int, uint64_t>;
	using View = Container::ContainerView;
	class NodeSliceIterator : std::iterator<std::forward_iterator_tag, std::pair<size_t, View>>
	{
		using map_iterator = typename std::unordered_map<size_t, MapItem>::iterator;
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
				return std::make_pair(nodeindex, slice->slices.getView(std::get<0>(info), std::get<1>(info), std::get<2>(info)));
			}
			else
			{
				auto nodeindex = mappos->first;
				auto start = std::get<0>(mappos->second);
				auto end = std::get<1>(mappos->second);
				return std::make_pair(nodeindex, slice->slices.getView(start, end, std::get<2>(mappos->second)));
			}
		}
		const std::pair<size_t, const View> operator*() const
		{
			if (indexPos != -1)
			{
				auto nodeindex = slice->activeVectorMapIndices[indexPos];
				auto info = (*(slice->vectorMap))[nodeindex];
				return std::make_pair(nodeindex, slice->slices.getView(std::get<0>(info), std::get<1>(info), std::get<2>(info)));
			}
			else
			{
				auto nodeindex = mappos->first;
				auto start = std::get<0>(mappos->second);
				auto end = std::get<1>(mappos->second);
				return std::make_pair(nodeindex, slice->slices.getView(start, end, std::get<2>(mappos->second)));
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
		using map_iterator = typename std::unordered_map<size_t, MapItem>::const_iterator;
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
				return std::make_pair(nodeindex, slice->slices.getView(std::get<0>(info), std::get<1>(info), std::get<2>(info)));
			}
			else
			{
				auto nodeindex = mappos->first;
				auto start = std::get<0>(mappos->second);
				auto end = std::get<1>(mappos->second);
				return std::make_pair(nodeindex, slice->slices.getView(start, end, std::get<2>(mappos->second)));
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
	vectorMap(nullptr)
	{
	}
	NodeSlice(std::vector<MapItem>* vectorMap) :
	vectorMap(vectorMap)
	{
	}
	void clearVectorMap()
	{
		assert(vectorMap != nullptr);
		for (auto index : activeVectorMapIndices)
		{
			(*vectorMap)[index] = std::make_tuple(0, 0, 0, 0);
		}
		activeVectorMapIndices.clear();
	}
	void reserve(size_t size)
	{
		slices.reserve(size);
	}
	void addNode(size_t nodeIndex, size_t size)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert(std::get<0>((*vectorMap)[nodeIndex]) == std::get<1>((*vectorMap)[nodeIndex]));
			(*vectorMap)[nodeIndex] = { slices.size(), slices.size() + size, 0, size };
			activeVectorMapIndices.push_back(nodeIndex);
		}
		else
		{
			assert(nodes.find(nodeIndex) == nodes.end());
			nodes[nodeIndex] = { slices.size(), slices.size() + size, 0, size };
		}
		slices.resize(slices.size() + size);
	}
	View node(size_t nodeIndex)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert(std::get<0>((*vectorMap)[nodeIndex]) != std::get<1>((*vectorMap)[nodeIndex]));
			return slices.getView(std::get<0>((*vectorMap)[nodeIndex]), std::get<1>((*vectorMap)[nodeIndex]), std::get<2>((*vectorMap)[nodeIndex]));
		}
		else
		{
			auto found = nodes.find(nodeIndex);
			assert(found != nodes.end());
			return slices.getView(std::get<0>(found->second), std::get<1>(found->second), std::get<2>(found->second));
		}
	}
	const View node(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert(std::get<0>((*vectorMap)[nodeIndex]) != std::get<1>((*vectorMap)[nodeIndex]));
			return slices.getView(std::get<0>((*vectorMap)[nodeIndex]), std::get<1>((*vectorMap)[nodeIndex]), std::get<2>((*vectorMap)[nodeIndex]));
		}
		else
		{
			auto found = nodes.find(nodeIndex);
			assert(found != nodes.end());
			return slices.getView(std::get<0>(found->second), std::get<1>(found->second), std::get<2>(found->second));
		}
	}
	bool hasNode(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			return std::get<0>((*vectorMap)[nodeIndex]) != std::get<1>((*vectorMap)[nodeIndex]);
		}
		else
		{
			return nodes.find(nodeIndex) != nodes.end();
		}
	}
	int minScore(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert(std::get<0>((*vectorMap)[nodeIndex]) != std::get<1>((*vectorMap)[nodeIndex]));
			return std::get<2>((*vectorMap)[nodeIndex]);
		}
		else
		{
			auto found = nodes.find(nodeIndex);
			assert(found != nodes.end());
			return std::get<2>(found->second);
		}
	}
	void setMinScoreIfSmaller(size_t nodeIndex, int score)
	{
		auto oldScore = minScore(nodeIndex);
		if (score < oldScore) setMinScore(nodeIndex, score);
	}
	void setMinScore(size_t nodeIndex, int score)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			std::get<2>((*vectorMap)[nodeIndex]) = score;
		}
		else
		{
			std::get<2>(nodes[nodeIndex]) = score;
		}
	}
	size_t startIndex(size_t nodeIndex) const
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			assert(std::get<0>((*vectorMap)[nodeIndex]) != std::get<1>((*vectorMap)[nodeIndex]));
			return std::get<3>((*vectorMap)[nodeIndex]);
		}
		else
		{
			auto found = nodes.find(nodeIndex);
			assert(found != nodes.end());
			return std::get<3>(found->second);
		}
	}
	void setStartIndex(size_t nodeIndex, size_t index)
	{
		if (vectorMap != nullptr)
		{
			assert(nodeIndex < vectorMap->size());
			std::get<3>((*vectorMap)[nodeIndex]) = index;
		}
		else
		{
			std::get<3>(nodes[nodeIndex]) = index;
		}
	}
	size_t size() const
	{
		if (vectorMap != nullptr)
		{
			return activeVectorMapIndices.size();
		}
		else
		{
			return nodes.size();
		}
	}
	NodeSliceIterator begin()
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceIterator { this, nodes.end(), 0 };
		}
		else
		{
			return NodeSliceIterator { this, nodes.begin() };
		}
	}
	NodeSliceIterator end()
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceIterator { this, nodes.end(), activeVectorMapIndices.size() };
		}
		else
		{
			return NodeSliceIterator { this, nodes.end() };
		}
	}
	NodeSliceConstIterator begin() const
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceConstIterator { this, nodes.end(), 0 };
		}
		else
		{
			return NodeSliceConstIterator { this, nodes.begin() };
		}
	}
	NodeSliceConstIterator end() const
	{
		if (vectorMap != nullptr)
		{
			return NodeSliceConstIterator { this, nodes.end(), activeVectorMapIndices.size() };
		}
		else
		{
			return NodeSliceConstIterator { this, nodes.end() };
		}
	}
	NodeSlice getFrozenSqrtEndScores() const
	{
		NodeSlice result;
		result.slices = slices.getFrozenSqrtEndScores();
		if (vectorMap != nullptr)
		{
			for (auto index : activeVectorMapIndices)
			{
				result.nodes[index] = (*vectorMap)[index];
			}
		}
		else
		{
			result.nodes = nodes;
		}
		return result;
	}
	NodeSlice getFrozenScores() const
	{
		NodeSlice result;
		result.slices = slices.getFrozenScores();
		if (vectorMap != nullptr)
		{
			for (auto index : activeVectorMapIndices)
			{
				result.nodes[index] = (*vectorMap)[index];
			}
		}
		else
		{
			result.nodes = nodes;
		}
		return result;
	}
private:
	std::vector<MapItem>* vectorMap;
	std::vector<size_t> activeVectorMapIndices;
	std::unordered_map<size_t, MapItem> nodes;
	Container slices;
	friend class NodeSliceIterator;
	friend class NodeSliceConstIterator;
};

#endif