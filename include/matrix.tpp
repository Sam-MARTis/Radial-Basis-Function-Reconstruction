#pragma once




template <typename T>
void MatrixSolver<T>::jacobiSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, const uint maxIterations, const T tolerance)
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

template <typename T>
void MatrixSolver<T>::gaussSeidelSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, const uint maxIterations, const T tolerance)
{
    uint n = A.numRows();
    assert(n == A.numColumns());
    bool converged = false;
    for (uint iter = 0; iter < maxIterations; ++iter)
    {
        std::vector<T> x_old = x;
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
            x[i] = (b[i] - sum) / A(i, i);
        }
        
        T norm = 0;
        for (uint i = 0; i < n; ++i)
        {
            norm += (x[i] - x_old[i]) * (x[i] - x_old[i]);
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


template <typename T>
bool MatrixSolver<T>::checkDiagonalDominance(const Matrix<T> &A)
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