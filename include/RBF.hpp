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
    Matrix<T> coefficientsDerivatives;
    std::vector<T> coefficients;
    public:
    RBFs(uint _numCells, uint _numFunctions, std::vector<Vec<T, N>>& _kernelCenters, std::unique_ptr<Function<T, N>> _function):
    numRBFs(_numFunctions), numCells(_numCells),
    kernelCenters(_kernelCenters), ConnectivityMatrix(_numCells, _numFunctions),
    coefficientsDerivatives(_numCells, _numFunctions), kernelFunction(std::move(_function)),
    coefficients(numRBFs)
    {
        /*
         *IMPORTANT
         *The number of rbs must be the same as the number of cells.
         *The indices must align with corresponding cells
         *The kth rbf corresponds to the kth cell!
         *Here even the rbf centers and cell centers align but that need not be the case.
         *In the future we can move kernel centers for some upwind like interpolation
         *BUT THE INDICES WILL STILL ALIGN. DONT FORGET THIS.
         */
        assert(_numCells == _numFunctions);
    }
    [[deprecated]] void computeConnectivityMatrix(const Mesh<T, N> &mesh);
    void computeLocalCellConnectivityMatrix(const Mesh<T, N> &mesh);
    void computeCoefficientDerivatives(uint maxIterations, T tolerance);
    [[deprecated]] void computeRBFCoefficients(const std::vector<T>& cellAveragesX, uint maxIterations, T tolerance);
    [[deprecated]]void computeRBFCoefficientsViaDerivatives(const std::vector<T>& cellAveragesX);

    void updateMeshRBFCoefficients(const Mesh<T, N> &mesh);
    T value(const Vec<T, N>& point) const;
    T valueLocal(const Vec<T, N>& point, const std::vector<Index>& activeRBFIndices, const std::vector<T>& coefficients) const;
    T integrate(T lb, T ub) const;
};

#include "RBF.tpp"
// class MultiQuadratic