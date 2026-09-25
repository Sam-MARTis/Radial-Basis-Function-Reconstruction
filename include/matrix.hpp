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
    static void jacobiSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, uint maxIterations, T tolerance)
    {
        uint n = A.numRows();
        assert(n == A.numColumns());
        assert(b.size() == n);
        assert(x.size() == n);
        std::vector<T> x_new(n, 0);
        bool converged = false;
        for (uint iter = 0; iter < maxIterations; ++iter)
        {
            for (uint i = 0; i < n; ++i)
            {
                T sum = 0;
                for (uint j = 0; j < n; ++j)
                {
                    if (j != i)
                    {
                        sum += A(i, j) * x[j];
                    }
                }
                x_new[i] = (b[i] - sum) / A(i, i);
            }

            // Check for convergence
            T norm = 0;
            for (uint i = 0; i < n; ++i)
            {
                norm += (x_new[i] - x[i]) * (x_new[i] - x[i]);
            }
            norm = std::sqrt(norm);

            if (norm < tolerance)
            {
                x = x_new;
                converged = true;
                break;
            }

            x = x_new;
        }
        if (!converged)
        {

            std::cerr << "Jacobi method did not converge within the maximum number of iterations." << std::endl;
        }
    }
    static void gaussSeidelSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, uint maxIterations, T tolerance)
    {
        const T overrelaxationFactor = 1.8;
        uint n = A.numRows();
        assert(n == A.numColumns());
        bool converged = false;
        for (uint iter = 0; iter < maxIterations; ++iter)
        {
            T norm = 0;
            for (uint i = 0; i < n; ++i)
            {
                T sum = 0;
                for (uint j = 0; j < n; ++j)
                {
                    if (j != i)
                    {
                        sum += A(i, j) * x[j];
                    }
                }
                T x_new = (b[i] - sum) / A(i, i);
                norm += (x_new - x[i]) * (x_new - x[i]);
                x[i] += overrelaxationFactor * (x_new - x[i]);
            }
            norm = std::sqrt(norm);

            if (norm < tolerance)
            {
                converged = true;
                break;
            }
        }
        if (!converged)
        {
            std::cerr << "Gauss-Seidel method did not converge within the maximum number of iterations." << std::endl;
        }
    }
    [[nodiscard]] static bool checkDiagonalDominance(const Matrix<T> &A)
    {
        uint n = A.numRows();
        assert(n == A.numColumns());
        for (uint i = 0; i < n; ++i)
        {
            T diagonalElement = std::abs(A(i, i));
            T sumOfOtherElements = 0;
            for (uint j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    sumOfOtherElements += std::abs(A(i, j));
                }
            }
            if (diagonalElement < sumOfOtherElements)
            {
                return false;
            }
        }
        return true;
    }
};


#include "matrix.tpp"