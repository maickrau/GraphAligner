#ifndef OrderedIndexKeeper_h
#define OrderedIndexKeeper_h

template <typename T>
class OrderedIndexKeeper
{
public:
	OrderedIndexKeeper()
	{

	}

	OrderedIndexKeeper(const std::vector<T>& initial) :
	zeroscore(std::numeric_limits<T>::max()),
	scores(),
	order(),
	reverseOrder(),
	breakpoints()
	{
		Initialize(initial);
	}

	void Initialize(const std::vector<T>& initial)
	{
		zeroscore = std::numeric_limits<T>::max();
		order.clear();
		scores.clear();
		reverseOrder.clear();
		breakpoints.clear();
		order.reserve(initial.size());
		scores.reserve(initial.size());
		reverseOrder.reserve(initial.size());
		breakpoints.reserve(initial.size()+4);
		for (size_t i = 0; i < initial.size(); i++)
		{
			zeroscore = std::min(zeroscore, initial[i] - 2);
		}
		//counting sort
		std::vector<T> counts;
		counts.resize(initial.size(), 0);
		for (size_t i = 0; i < initial.size(); i++)
		{
			assert(initial[i] >= (zeroscore + 2));
			assert(initial[i] < (zeroscore + 2) + initial.size());
			counts[initial[i] - (zeroscore + 2)]++;
		}
		for (size_t i = 1; i < initial.size(); i++)
		{
			counts[i] += counts[i-1];
		}
		order.resize(initial.size(), std::numeric_limits<T>::max());
		for (size_t i = 0; i < initial.size(); i++)
		{
			assert(counts[initial[i]-(zeroscore + 2)] > 0);
			assert(counts[initial[i]-(zeroscore + 2)] - 1 < initial.size());
			assert(order[counts[initial[i]-(zeroscore + 2)] - 1] == std::numeric_limits<T>::max());
			order[counts[initial[i]-(zeroscore + 2)] - 1] = i;
			counts[initial[i]-(zeroscore + 2)]--;
		}
#ifndef NDEBUG
		for (size_t i = 0; i < initial.size(); i++)
		{
			assert(order[i] != std::numeric_limits<T>::max());
		}
#endif
		breakpoints.push_back(0);
		breakpoints.push_back(0);
		breakpoints.push_back(0);
		scores.resize(initial.size(), std::numeric_limits<T>::max());
		reverseOrder.resize(initial.size(), std::numeric_limits<T>::max());
		for (size_t i = 0; i < order.size(); i++)
		{
			assert(initial[order[i]] >= (zeroscore + 2));
			assert(initial[order[i]] < (zeroscore + 2) + initial.size());
			scores[i] = initial[order[i]];
			if (scores[i] > scores[breakpoints.back()])
			{
				breakpoints.push_back(i);
			}
			reverseOrder[order[i]] = i;
		}
		assert(breakpoints.size() <= initial.size()+2);
		breakpoints.resize(initial.size()+4, initial.size());
		assert(order.size() == scores.size());
		assert(scores.size() == reverseOrder.size());
		assert(breakpoints.size() == reverseOrder.size()+4);
#ifndef NDEBUG
		for (size_t i = 0; i < initial.size(); i++)
		{
			assert(scores[i] != std::numeric_limits<T>::max());
			assert(reverseOrder[i] != std::numeric_limits<T>::max());
			assert(order[i] != std::numeric_limits<T>::max());
			assert(reverseOrder[order[i]] == i);
			assert(order[reverseOrder[i]] == i);
		}
#endif
	}

	size_t size() const
	{
		return order.size();
	}

	T GetColumn(T position) const
	{
		return order[position];
	}

	T GetPosition(T column) const
	{
		return reverseOrder[column];
	}

	T GetValue(T column) const
	{
		return scores[GetPosition(column)];
	}

	void SetIfSmaller(T column, T value)
	{
		auto position = GetPosition(column);
		auto oldValue = scores[position];
		if (value < oldValue)
		{
			for (T i = 0; i < oldValue - value; i++)
			{
				decrement(position);
				position = GetPosition(column);
				assert(GetColumn(position) == column);
			}
		}
		assert(GetValue(column) <= value);
	}

	void SetValue(T column, T value)
	{
		auto position = GetPosition(column);
		auto oldValue = scores[position];
		if (value < oldValue)
		{
			for (T i = 0; i < oldValue - value; i++)
			{
				decrement(position);
				position = GetPosition(column);
				assert(GetColumn(position) == column);
			}
		}
		if (value > oldValue)
		{
			for (T i = 0; i < value - oldValue; i++)
			{
				increment(position);
				position = GetPosition(column);
				assert(GetColumn(position) == column);
			}
		}
		assert(GetValue(column) == value);
	}

	void FixZeroScore()
	{
		auto oldZeroScore = zeroscore;
		auto newZeroScore = scores[0]-2;
		if (newZeroScore == oldZeroScore) return;
		zeroscore = newZeroScore;
		breakpoints[0] = 0;
		breakpoints[1] = 0;
		breakpoints[2] = 0;
		size_t breakpointpos = 3;
		for (size_t i = 1; i < scores.size(); i++)
		{
			if (scores[i] > scores[i-1])
			{
				for (; zeroscore + breakpointpos <= scores[i]; breakpointpos++)
				{
					assert(breakpointpos < breakpoints.size());
					breakpoints[breakpointpos] = i;
				}
			}
		}
		for (; breakpointpos < breakpoints.size(); breakpointpos++)
		{
			breakpoints[breakpointpos] = scores.size();
		}
	}

private:

	void increment(T position)
	{
		assert(position >= 0);
		assert(position < scores.size());
		auto score = scores[position];
		assert(score >= zeroscore);
		assert(score - zeroscore + 1 < breakpoints.size());
		auto blockEnd = breakpoints[score-zeroscore+1]-1;
		assert(blockEnd >= 0);
		assert(blockEnd < order.size());
		auto oldReverseOrder = order[position];
		auto endReverseOrder = order[blockEnd];
		std::swap(order[position], order[blockEnd]);
		scores[blockEnd]++;
		breakpoints[score-zeroscore+1]--;
		std::swap(reverseOrder[oldReverseOrder], reverseOrder[endReverseOrder]);
	}

	void decrement(T position)
	{
		assert(position >= 0);
		assert(position < scores.size());
		auto score = scores[position];
		assert(score >= zeroscore);
		assert(score - zeroscore + 1 < breakpoints.size());
		auto blockStart = breakpoints[score-zeroscore];
		assert(blockStart >= 0);
		assert(blockStart < order.size());
		auto oldReverseOrder = order[position];
		auto startReverseOrder = order[blockStart];
		std::swap(order[position], order[blockStart]);
		scores[blockStart]--;
		breakpoints[score-zeroscore]++;
		std::swap(reverseOrder[oldReverseOrder], reverseOrder[startReverseOrder]);
	}

	T zeroscore;
	std::vector<T> scores; //S
	std::vector<T> order; //O
	std::vector<T> reverseOrder; //R
	std::vector<T> breakpoints; //B
};

#endif