#ifndef Array2D_h
#define Array2D_h

#include <cassert>

template<bool columnMajor>
size_t getIndex(size_t row, size_t column, size_t rows, size_t columns)
{
	assert(false);
	return 0;
}

template<>
size_t getIndex<false>(size_t row, size_t column, size_t rows, size_t columns)
{
	return row*columns+column;
}

template<>
size_t getIndex<true>(size_t row, size_t column, size_t rows, size_t columns)
{
	return column*rows+row;
}

template <typename T, bool columnMajor>
class Array2D
{
public:
	Array2D(size_t rows, size_t columns, T value) :
	rows(rows),
	columns(columns),
	data(nullptr)
	{
		data = new T[rows*columns];
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < columns; j++)
			{
				get(i, j) = value;
			}
		}
	};
	Array2D(const Array2D& other) :
	rows(other.rows),
	columns(other.columns),
	data(nullptr)
	{
		copy(other);
	};
	Array2D(Array2D&& other) :
	rows(other.rows),
	columns(other.columns),
	data(other.data)
	{
		other.data = nullptr;
	};
	~Array2D()
	{
		deallocate();
	}
	Array2D& operator=(const Array2D& other)
	{
		if (&other == this) return *this;
		deallocate();
		copy(other);
		return *this;
	};
	Array2D& operator=(Array2D&& other)
	{
		if (&other == this) return *this;
		deallocate();
		data = other.data;
		other.data = nullptr;
		return *this;
	};
	T& operator()(size_t row, size_t column)
	{
		return get(row, column);
	};
	T operator()(size_t row, size_t column) const
	{
		return get(row, column);
	}
	T& get(size_t row, size_t column)
	{
		return data[getIndex<columnMajor>(row, column, rows, columns)];
	}
	T get(size_t row, size_t column) const
	{
		return data[getIndex<columnMajor>(row, column, rows, columns)];
	}
	size_t sizeColumns() const
	{
		return columns;
	}
	size_t sizeRows() const
	{
		return rows;
	}
private:
	// size_t getIndex(size_t row, size_t column) const
	// {
	// 	assert(false);
	// 	return 0;
	// }
	void deallocate()
	{
		if (data != nullptr) delete [] data;
		data = nullptr;
	}
	void copy(const Array2D& other)
	{
		deallocate();
		data = new T[rows*columns];
		memcpy(data, other.data, sizeof(T)*rows*columns);
	}
	size_t rows;
	size_t columns;
	T* data;
};

#endif