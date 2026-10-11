#pragma once

#include <vector>
#include "vec.hpp"
#include "constants.hpp"
#include "matrix.hpp"


template <typename T, uint N>
struct EdgeProperty
{
    // std::array<T, static_cast<uint>(pow(2, N-1))> vertices;

    std::array<T, N> normal;
    Vec<T, N+2> numericalFlux;
    std::pair<Vec<T, N+2>, Vec<T, N+2>> reconstructedPushedU; // If normal . direction > 0, use element 1, else 0.
    std::pair<Index, Index> connectingCells;
    T measure = 1;
};


template <typename T, uint N>
struct Edges
{
    std::vector<EdgeProperty<T, N>> edgeProperty;
    std::vector<std::array<Vec<T, N>, static_cast<uint>(1) << (N - 1)>> vertices;
    std::vector<Vec<T, N>> edgeEvaluationPoints;
};

template <typename T, uint N>
struct Neighbourhood
{
    /*
     * Two neighbourhoods per cell, with independent depths:
     *
     *  - neighbourCellIndices / ALocalInv / rbfCoefficients size the local RBF fit.
     *  - extremaNeighbourCellIndices feeds the min/max used by the limiter bounds.
     *
     * They are separate because a TVD limiter wants the immediate stencil while the
     * RBF fit needs a wider support. Both include the cell's own index.
     */
    Index numNeighbours = 0;
    std::vector<Index> neighbourCellIndices;
    std::array<std::vector<T>, N+2> rbfCoefficients;
    Matrix<T> ALocalInv;

    Index numExtremaNeighbours = 0;
    std::vector<Index> extremaNeighbourCellIndices;
    std::pair<Vec<T, N+2>, Vec<T, N+2>> extremeValues;
};

template <typename T, uint N>
struct Cells
{
    std::vector<std::vector<Index>> edgeIndices;
    std::vector<Vec<T, N+2>> cellAverages;
    // std::vector<Vec<Vec<T, N>, N+2>> reconstructedCellGradients;
    std::vector<Neighbourhood<T, N>> neighbourhoods;
    std::vector<T> cellVolumes;
    std::vector<Vec<T, N>> cellCentersPositions;
};

