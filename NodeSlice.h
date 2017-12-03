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
			WordContainerIterator result;
			result.container = container;
			result.pos = pos + off;
			return result;
		}
		WordContainerIterator operator-(size_t off) const
		{
			WordContainerIterator result;
			result.container = container;
			result.pos = pos - off;
			return result;
		}
		WordContainerIterator& operator+=(size_t off)
		{
			pos += off;
		}
		WordContainerIterator& operator-=(size_t off)
		{
			pos -= off;
		}
		size_t operator-(const WordContainerIterator& other) const
		{
			return pos - other.pos;
		}
	private:
		WordContainer* container;
		size_t pos;
	};
	using iterator = WordContainerIterator;
	using const_iterator = const WordContainerIterator;
	WordContainer() :
	minEndScore(0),
	minStartScore(0),
	frozen(0)
	{
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
	void freezeScores()
	{
		assert(frozen == 0);
		frozenSlices.resize(mutableSlices.size());
		minStartScore = mutableSlices[0].scoreBeforeStart;
		for (size_t i = 1; i < mutableSlices.size(); i++)
		{
			minStartScore = std::min(minStartScore, mutableSlices[i].scoreBeforeStart);
		}
		for (size_t i = 0; i < mutableSlices.size(); i++)
		{
			frozenSlices[i].VP = mutableSlices[i].VP;
			frozenSlices[i].VN = mutableSlices[i].VN;
			assert(mutableSlices[i].scoreBeforeStart >= minStartScore);
			assert(mutableSlices[i].scoreBeforeStart - minStartScore < std::numeric_limits<decltype(frozenSlices[i].plusMinScore)>::max());
			frozenSlices[i].plusMinScore = mutableSlices[i].scoreBeforeStart - minStartScore;
		}
		{
			std::vector<WordSlice> empty;
			std::swap(mutableSlices, empty);
		}
		frozen = 1;
	}
	void freezeSqrtEndscores()
	{
		assert(frozen == 0);
		frozenSqrtSlices.resize(mutableSlices.size());
		minEndScore = mutableSlices[0].scoreEnd;
		for (size_t i = 1; i < mutableSlices.size(); i++)
		{
			minEndScore = std::min(minEndScore, mutableSlices[i].scoreEnd);
		}
		for (size_t i = 0; i < mutableSlices.size(); i++)
		{
			frozenSqrtSlices[i].VPVNLastBit = 0;
			frozenSqrtSlices[i].VPVNLastBit |= mutableSlices[i].VP >> 63;
			frozenSqrtSlices[i].VPVNLastBit |= (mutableSlices[i].VN >> 62) & 2;
			frozenSqrtSlices[i].VPVNLastBit |= mutableSlices[i].scoreEndExists << 2;
			assert(mutableSlices[i].scoreEnd >= minEndScore);
			assert(mutableSlices[i].scoreEnd - minEndScore < std::numeric_limits<decltype(frozenSqrtSlices[i].plusMinScore)>::max());
			frozenSqrtSlices[i].plusMinScore = mutableSlices[i].scoreEnd - minEndScore;
		}
		{
			std::vector<WordSlice> empty;
			std::swap(mutableSlices, empty);
		}
		frozen = 2;
	}
	iterator begin()
	{
		return iterator { this, 0 };
	}
	iterator end()
	{
		return iterator { this, size() };
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
	void addNode(size_t nodeIndex, size_t size)
	{
		nodes[nodeIndex].resize(size);
	}
	Container& node(size_t nodeIndex)
	{
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return found->second;
	}
	const Container& node(size_t nodeIndex) const
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
		auto found = nodes.find(nodeIndex);
		assert(found != nodes.end());
		return found->second.minScore;
	}
	void setMinScore(size_t nodeIndex, int score)
	{
		nodes[nodeIndex].minScore = score;
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
	void freezeSqrtEndScores()
	{
		for (auto& pair : nodes)
		{
			pair.second.freezeSqrtEndscores();
		}
	}
	void freezeScores()
	{
		for (auto& pair : nodes)
		{
			pair.second.freezeScores();
		}
	}
private:
	std::unordered_map<size_t, Container> nodes;
};

#endif