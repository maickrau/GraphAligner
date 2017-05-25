#include <cstdlib>
#include "SparseBoolMatrix.h"

SparseBoolMatrix::SparseBoolMatrix(size_t numColumns, size_t numRows) :
numColumns(numColumns),
numRows(numRows)
{
	rows.resize(numRows);
}

bool SparseBoolMatrix::get(size_t column, size_t row) const
{
	return rows[row].count(column) > 0;
}

void SparseBoolMatrix::set(size_t column, size_t row)
{
	rows[row].insert(column);
}

void SparseBoolMatrix::unset(size_t column, size_t row)
{
	rows[row].erase(column);
}

bool SparseBoolMatrix::operator()(size_t column, size_t row) const
{
	return get(column, row);
}

std::set<size_t>::const_iterator SparseBoolMatrix::rowStart(size_t row) const
{
	return rows[row].begin();
}

std::set<size_t>::const_iterator SparseBoolMatrix::rowEnd(size_t row) const
{
	return rows[row].end();
}

size_t SparseBoolMatrix::sizeRows() const
{
	return numRows;
}

size_t SparseBoolMatrix::sizeColumns() const
{
	return numColumns;
}

size_t SparseBoolMatrix::rowSize(size_t row) const
{
	return rows[row].size();
}