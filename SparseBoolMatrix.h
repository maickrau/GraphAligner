#ifndef SparseBoolMatrix_h
#define SparseBoolMatrix_h

#include <unordered_set>
#include <vector>

class SparseBoolMatrix
{
public:
	SparseBoolMatrix(size_t columns, size_t rows);
	bool get(size_t column, size_t row) const;
	void set(size_t column, size_t row);
	void unset(size_t column, size_t row);
	bool operator()(size_t column, size_t row) const;
	std::unordered_set<size_t>::const_iterator rowStart(size_t row) const;
	std::unordered_set<size_t>::const_iterator rowEnd(size_t row) const;
	template<typename Iterator>
	void addRow(size_t row, Iterator start, Iterator end)
	{
		rows[row].insert(start, end);
	}
	size_t sizeRows() const;
	size_t sizeColumns() const;
	size_t rowSize(size_t row) const;
private:
	size_t numColumns;
	size_t numRows;
	std::vector<std::unordered_set<size_t>> rows;
};

#endif
