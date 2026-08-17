#pragma once




template <typename T>
void MatrixSolver<T>::jacobiSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, const uint maxIterations, const T tolerance)
{
    uint n = A.numRows();
    assert(n == A.numColumns());
    std::vector<T> x_new(n, 0);
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
            break;
        }

        x = x_new;
    }
}

template <typename T>
void MatrixSolver<T>::gaussSeidelSolver(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &x, const uint maxIterations, const T tolerance)
{
    uint n = A.numRows();
    assert(n == A.numColumns());
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
            break;
        }
    }
}