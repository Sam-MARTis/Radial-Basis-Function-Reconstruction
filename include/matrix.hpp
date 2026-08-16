#pragma once
#include "constants.hpp"
#include <vector>

template <typename T>
class Matrix
{
    uint rows;
    uint columns;
    std::vector<T> data;
public:
    Matrix(const uint _rows, const uint _columns): rows(_rows), columns(_columns), data(rows * columns) {}
    T& operator()(const uint row, const uint column)
    {
        return data[row * columns + column];
    }
    const T& operator()(const uint row, const uint column) const
    {
        return data[row * columns + column];
    }
    [[nodiscard]] uint numRows() const { return rows; }
    [[nodiscard]] uint numColumns() const { return columns; }
};