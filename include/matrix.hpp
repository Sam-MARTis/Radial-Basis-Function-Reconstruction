#pragma once
#include "constants.hpp"
#include <vector>
#include <iostream>

template <typename T>
class Matrix
{
    uint rows;
    uint columns;
public:
    std::vector<T> data;
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

template <typename T>
class MatrixSolver
{
public:
    static void jacobiSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, uint maxIterations, T tolerance);
    static void gaussSeidelSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, uint maxIterations, T tolerance);
    [[nodiscard]] static bool checkDiagonalDominance(const Matrix<T> &A);
};


#include "matrix.tpp"