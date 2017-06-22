#ifndef SparseBoolMatrix_h
#define SparseBoolMatrix_h

#include <vector>
#include "ThreadReadAssertion.h"

template <typename ContainerType>
class SparseBoolMatrix
{
public:
	SparseBoolMatrix(size_t numColumns, size_t numRows) :
	numColumns(numColumns),
	numRows(numRows)
	{
		rows.resize(numRows);
	}
	size_t getSolidIndex(size_t column, size_t row) const
	{
		return rows[row].getSolidIndex(column);
	}
	bool get(size_t column, size_t row) const
	{
		return rows[row].count(column) > 0;
	}
	void set(size_t column, size_t row)
	{
		rows[row].insert(column);
	}
	void setBlock(size_t columnStart, size_t columnEnd, size_t row)
	{
		rows[row].insertBlock(columnStart, columnEnd);
	}
	void unset(size_t column, size_t row)
	{
		rows[row].erase(column);
	}
	bool operator()(size_t column, size_t row) const
	{
		return get(column, row);
	}
	typename ContainerType::const_iterator rowStart(size_t row) const
	{
		return rows[row].begin();
	}
	typename ContainerType::const_iterator rowEnd(size_t row) const
	{
		return rows[row].end();
	}
	template<typename Iterator>
	void addRow(size_t row, Iterator start, Iterator end)
	{
		rows[row].insert(start, end);
	}
	size_t sizeRows() const
	{
		return numRows;
	}
	size_t sizeColumns() const
	{
		return numColumns;
	}
	size_t rowSize(size_t row) const
	{
		return rows[row].size();
	}
	size_t rowSizeBlock(size_t row) const
	{
		return rows[row].numBlocks();
	}
	size_t totalOnes() const
	{
		size_t result = 0;
		for (size_t i = 0; i < rows.size(); i++)
		{
			result += rows[i].size();
		}
		return result;
	}
private:
	size_t numColumns;
	size_t numRows;
	std::vector<ContainerType> rows;
};

#endif
