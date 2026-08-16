#pragma once
#include <memory>

#include "function.hpp"
#include "matrix.hpp"





template<typename T, uint N>
class RBFs{
    uint numRBFs = 0;
    std::unique_ptr<Function<T, N>> kernelFunction;
    std::vector<Vec<T, N>> kernelCenters;
    Matrix<T> ConnectivityMatrix;
    public:
    RBFs(uint numCells, uint numFunctions, std::unique_ptr<Function<T, N>> function): numRBFs(numFunctions), ConnectivityMatrix(numFunctions, numCells), kernelFunction(std::move(function)) {}
    void computeConnectivity();
};


// class MultiQuadratic