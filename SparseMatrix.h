#ifndef SparseMatrix_h
#define SparseMatrix_h

#include <unordered_map>
#include <vector>
#include "ThreadReadAssertion.h"

template <typename T, typename BandMatrix>
class SparseMatrix
{
public:
	SparseMatrix(size_t numColumns, size_t numRows, const BandMatrix& band) :
	numColumns(numColumns),
	numRows(numRows),
	band(band)
	{
		rows.resize(numRows);
	}
	void set(size_t column, size_t row, T value)
	{
		assert(row < rows.size());
		if (rows[row].size() == 0) rows[row].resize(band.rowSize(row));
		auto solidIndex = band.getSolidIndex(column, row);
		assert(solidIndex < rows[row].size());
		rows[row][solidIndex] = value;
	}
	T operator()(size_t column, size_t row) const
	{
		return get(column, row);
	}
	T get(size_t column, size_t row) const
	{
		assert(row < rows.size());
		auto solidIndex = band.getSolidIndex(column, row);
		assert(solidIndex < rows[row].size());
		return rows[row][solidIndex];
	}
	bool exists(size_t column, size_t row) const
	{
		assert(row < rows.size());
		//todo fix
		//this will break on an assertion if it does not exist, 
		//but this function is only called inside assertions anyway
		auto solidIndex = band.getSolidIndex(column, row);
		assert(solidIndex < rows[row].size());
		return solidIndex != -1;
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
	const BandMatrix& band;
	size_t numColumns;
	size_t numRows;
	std::vector<std::vector<T>> rows;
};

#endif