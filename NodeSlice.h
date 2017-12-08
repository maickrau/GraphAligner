#ifndef NodeSlice_h
#define NodeSlice_h

#include <limits>
#include <unordered_map>
#include <vector>
#include "ThreadReadAssertion.h"

template <typename LengthType, typename ScoreType, typename Word>
class WordContainer
{
public:
	class RowConfirmation
	{
	public:
		RowConfirmation(char rows, bool partial) : rows(rows), partial(partial) {};
		char rows;
		bool partial;
		bool operator>(const RowConfirmation& other) const
		{
			return rows > other.rows || (rows == other.rows && partial && !other.partial);
		}
		bool operator<(const RowConfirmation& other) const
		{
			return rows < other.rows || (rows == other.rows && !partial && other.partial);
		}
		bool operator==(const RowConfirmation& other) const
		{
			return rows == other.rows && partial == other.partial;
		}
		bool operator!=(const RowConfirmation& other) const
		{
			return !(*this == other);
		}
		bool operator>=(const RowConfirmation& other) const
		{
			return !(*this < other);
		}
		bool operator<=(const RowConfirmation& other) const
		{
			return !(*this > other);
		}
	};
	class WordSlice
	{
	public:
		WordSlice() :
		VP(0),
		VN(0),
		scoreEnd(0),
		scoreBeforeStart(0),
		confirmedRows(0, false),
		scoreBeforeExists(false),
		scoreEndExists(true)
		{}
		WordSlice(Word VP, Word VN, ScoreType scoreEnd, ScoreType scoreBeforeStart, int confirmedRows, bool scoreBeforeExists) :
		VP(VP),
		VN(VN),
		scoreEnd(scoreEnd),
		scoreBeforeStart(scoreBeforeStart),
		confirmedRows(confirmedRows, false),
		scoreBeforeExists(scoreBeforeExists),
		scoreEndExists(true)
		{}
		Word VP;
		Word VN;
		ScoreType scoreEnd;
		ScoreType scoreBeforeStart;
		RowConfirmation confirmedRows;
		bool scoreBeforeExists;
		bool scoreEndExists;
	};
	class SmallSlice
	{
	public:
		Word VP;
		Word VN;
		uint16_t plusMinScore;
	};
	class TinySlice
	{
	public:
		uint16_t plusMinScore;
		char VPVNLastBit;
	};
	class WordContainerIterator : std::iterator<std::random_access_iterator_tag, WordSlice>
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
		WordSlice& operator*()
		{
			return (*container)[pos];
		}
		const WordSlice operator*() const
		{
			return (*container)[pos];
		}
		WordSlice& operator[](size_t index)
		{
			return (*container)[pos + index];
		}
		const WordSlice operator[](size_t index) const
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
	class WordContainerConstIterator : std::iterator<std::random_access_iterator_tag, WordSlice>
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
		const WordSlice operator*() const
		{
			return (*container)[pos];
		}
		const WordSlice operator[](size_t index) const
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
		WordSlice& back()
		{
			return *(_end - 1);
		}
		const WordSlice back() const
		{
			return *(_cend - 1);
		}
		size_t size() const
		{
			return _cend - _cbegin;
		}
		WordSlice& operator[](size_t pos)
		{
			return *(_begin + pos);
		}
		const WordSlice operator[](size_t pos) const
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
	frozen(0)
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
	WordSlice& operator[](size_t index)
	{
		assert(frozen == 0);
		return mutableSlices[index];
	}
	const WordSlice operator[](size_t index) const
	{
		if (frozen == 1)
		{
			return { frozenSlices[index].VP, frozenSlices[index].VN, 0, minStartScore + frozenSlices[index].plusMinScore, 64, false };
		}
		else if (frozen == 2)
		{
			bool VP = frozenSqrtSlices[index].VPVNLastBit & 1;
			bool VN = frozenSqrtSlices[index].VPVNLastBit & 2;
			WordSlice result { (Word)VP << 63, (Word)VN << 63, minEndScore + frozenSqrtSlices[index].plusMinScore, minEndScore + frozenSqrtSlices[index].plusMinScore - (VP ? 1 : 0) + (VN ? 1 : 0), 64, false };
			result.scoreEndExists = frozenSqrtSlices[index].VPVNLastBit & 4;
			return result;
		}
		else
		{
			return mutableSlices[index];
		}
	}
	WordContainer getFrozenScores() const
	{
		if (frozen == 1) return *this;
		assert(frozen == 0);
		WordContainer result;
		result.frozen = 1;
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
		}
		return result;
	}
	WordContainer getFrozenSqrtEndScores() const
	{
		if (frozen == 2) return *this;
		assert(frozen == 0);
		WordContainer result;
		result.frozen = 2;
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
			result.frozenSqrtSlices[i].VPVNLastBit |= mutableSlices[i].scoreEndExists << 2;
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
		return frozen == 0 ? mutableSlices.size() : frozen == 1 ? frozenSlices.size() : frozenSqrtSlices.size();
	}
	void resize(size_t size)
	{
		assert(frozen == 0);
		mutableSlices.resize(size);
	}
	WordSlice& back()
	{
		assert(frozen == 0);
		return mutableSlices.back();
	}
	const WordSlice back() const
	{
		return (*this)[size()-1];
	}
	void reserve(size_t size)
	{
		assert(frozen == 0);
		mutableSlices.reserve(size);
	}
	ScoreType minScore;
private:
	ScoreType minEndScore;
	ScoreType minStartScore;
	int frozen;
	std::vector<WordSlice> mutableSlices;
	std::vector<SmallSlice> frozenSlices;
	std::vector<TinySlice> frozenSqrtSlices;
};

template <typename T>
class NodeSlice
{
public:
	using Container = WordContainer<size_t, int, uint64_t>;
	using View = Container::ContainerView;
	class NodeSliceIterator : std::iterator<std::forward_iterator_tag, std::pair<size_t, View>>
	{
		using map_iterator = typename std::unordered_map<size_t, std::tuple<size_t, size_t, int>>::iterator;
	public:
		NodeSliceIterator(NodeSlice* slice, map_iterator pos) :
		slice(slice),
		mappos(pos)
		{
		}
		std::pair<size_t, View> operator*()
		{
			auto nodeindex = mappos->first;
			auto start = std::get<0>(mappos->second);
			auto end = std::get<1>(mappos->second);
			return std::make_pair(nodeindex, slice->slices.getView(start, end, std::get<2>(mappos->second)));
		}
		const std::pair<size_t, View> operator*() const
		{
			auto nodeindex = mappos->first;
			auto start = std::get<0>(mappos->second);
			auto end = std::get<1>(mappos->second);
			return std::make_pair(nodeindex, slice->slices.getView(start, end, std::get<2>(mappos->second)));
		}
		NodeSliceIterator& operator++()
		{
			++mappos;
			return *this;
		}
		bool operator==(const NodeSliceIterator& other) const
		{
			return slice == other.slice && mappos == other.mappos;
		}
		bool operator!=(const NodeSliceIterator& other) const
		{
			return !(*this == other);
		}
	private:
		NodeSlice* slice;
		map_iterator mappos;
	};
	class NodeSliceConstIterator : std::iterator<std::forward_iterator_tag, const std::pair<size_t, View>>
	{
		using map_iterator = typename std::unordered_map<size_t, std::tuple<size_t, size_t, int>>::const_iterator;
	public:
		NodeSliceConstIterator(const NodeSlice* slice, map_iterator pos) :
		slice(slice),
		mappos(pos)
		{
		}
		const std::pair<size_t, View> operator*() const
		{
			auto nodeindex = mappos->first;
			auto start = std::get<0>(mappos->second);
			auto end = std::get<1>(mappos->second);
			return std::make_pair(nodeindex, slice->slices.getView(start, end, std::get<2>(mappos->second)));
		}
		NodeSliceConstIterator& operator++()
		{
			++mappos;
			return *this;
		}
		bool operator==(const NodeSliceConstIterator& other) const
		{
			return slice == other.slice && mappos == other.mappos;
		}
		bool operator!=(const NodeSliceConstIterator& other) const
		{
			return !(*this == other);
		}
	private:
		const NodeSlice* slice;
		map_iterator mappos;
	};
	void reserve(size_t size)
	{
		slices.reserve(size);
	}
	void addNode(size_t nodeIndex, size_t size)
	{
		nodes[nodeIndex] = { slices.size(), slices.size() + size, 0 };
		slices.resize(slices.size() + size);
	}
	View node(size_t nodeIndex)
	{
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return slices.getView(std::get<0>(found->second), std::get<1>(found->second), std::get<2>(found->second));
	}
	const View node(size_t nodeIndex) const
	{
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return slices.getView(std::get<0>(found->second), std::get<1>(found->second), std::get<2>(found->second));
	}
	bool hasNode(size_t nodeIndex) const
	{
		return nodes.find(nodeIndex) != nodes.end();
	}
	int minScore(size_t nodeIndex) const
	{
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return std::get<2>(found->second);
	}
	void setMinScore(size_t nodeIndex, int score)
	{
		std::get<2>(nodes[nodeIndex]) = score;
	}
	size_t size() const
	{
		return nodes.size();
	}
	NodeSliceIterator begin()
	{
		return NodeSliceIterator { this, nodes.begin() };
	}
	NodeSliceIterator end()
	{
		return NodeSliceIterator { this, nodes.end() };
	}
	NodeSliceConstIterator begin() const
	{
		return NodeSliceConstIterator { this, nodes.begin() };
	}
	NodeSliceConstIterator end() const
	{
		return NodeSliceConstIterator { this, nodes.end() };
	}
	NodeSlice getFrozenSqrtEndScores() const
	{
		NodeSlice result;
		result.slices = slices.getFrozenSqrtEndScores();
		result.nodes = nodes;
		return result;
	}
	NodeSlice getFrozenScores() const
	{
		NodeSlice result;
		result.slices = slices.getFrozenScores();
		result.nodes = nodes;
		return result;
	}
private:
	std::unordered_map<size_t, std::tuple<size_t, size_t, int>> nodes;
	Container slices;
};

#endif