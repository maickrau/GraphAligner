#ifndef SparseBoolMatrix_h
#define SparseBoolMatrix_h

#include <vector>

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
	bool get(size_t column, size_t row) const
	{
		return rows[row].count(column) > 0;
	}
	void set(size_t column, size_t row)
	{
		rows[row].insert(column);
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
private:
	size_t numColumns;
	size_t numRows;
	std::vector<ContainerType> rows;
};

#endif
