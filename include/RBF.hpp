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
    std::vector<Vec<T, N>> kernelCenters;
    Matrix<T> ConnectivityMatrix;
    std::vector<T> coefficients;
    public:
    RBFs(uint _numCells, uint _numFunctions, std::unique_ptr<Function<T, N>> _function): numRBFs(_numFunctions), numCells(_numCells), ConnectivityMatrix(_numFunctions, _numCells), kernelFunction(std::move(_function)), coefficients(numRBFs) {}
    void computeConnectivityMatrix(const Mesh<T, N> &mesh);
    void computeRBFCoefficients(const std::vector<T>& cellAveragesX);
    T value(const Vec<T, N>& point);
};


// class MultiQuadratic