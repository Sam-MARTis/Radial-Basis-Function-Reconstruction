#pragma once
#include <memory>

#include "function.hpp"
#include "matrix.hpp"
#include "mesh.hpp"




template<typename T, uint N>
class RBFs{
    const uint numRBFs;
    const uint numCells;
    std::unique_ptr<Function<T, N>> kernelFunction;
    std::vector<Vec<T, N>>& kernelCenters;
    Matrix<T> ConnectivityMatrix;
    std::vector<T> coefficients;
    public:
    RBFs(uint _numCells, uint _numFunctions, std::vector<Vec<T, N>>& _kernelCenters, std::unique_ptr<Function<T, N>> _function): numRBFs(_numFunctions), numCells(_numCells), kernelCenters(_kernelCenters), ConnectivityMatrix(_numCells, _numFunctions), kernelFunction(std::move(_function)), coefficients(numRBFs) {}
    void computeConnectivityMatrix(const Mesh<T, N> &mesh);
    void computeRBFCoefficients(const std::vector<T>& cellAveragesX, const uint maxIterations, const T tolerance);
    T value(const Vec<T, N>& point) const;
    T integrate(T lb, T ub) const;
};

#include "RBF.tpp"
// class MultiQuadratic