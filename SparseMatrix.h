#ifndef SparseMatrix_h
#define SparseMatrix_h

#include <unordered_map>
#include <vector>

template <typename T>
class SparseMatrix
{
public:
	SparseMatrix(size_t numColumns, size_t numRows) :
	numColumns(numColumns),
	numRows(numRows)
	{
		rows.resize(numRows);
	}
	void set(size_t column, size_t row, T value)
	{
		rows[row][column] = value;
	}
	T operator()(size_t column, size_t row) const
	{
		return get(column, row);
	}
	T get(size_t column, size_t row) const
	{
		return rows[row].at(column);
	}
	bool exists(size_t column, size_t row) const
	{
		return rows[row].count(column) > 0;
	}
	size_t sizeColumns() const
	{
		return numColumns;
	}
	size_t sizeRows() const
	{
		return numRows;
	}
private:
	size_t numColumns;
	size_t numRows;
	std::vector<std::unordered_map<size_t, T>> rows;
};

#endif