#pragma once

#include <vector>
#include "vec.hpp"
#include "constants.hpp"



template <typename T, uint N>
struct EdgeProperty
{
    // std::array<T, static_cast<uint>(pow(2, N-1))> vertices;

    std::array<T, N> normal;
    std::array<T, 2> normalVelFlux; // If normal . direction > 0, use element 1, else 0.
    T numericalFlux;
    std::pair<uint, uint> connectingCells;
};

template <typename T, uint N>
struct Edges
{
    std::vector<EdgeProperty<T, N>> edgeProperty;
    std::vector<std::array<Vec<T, N>, static_cast<uint>(1) << (N - 1)>> vertices;
};

template <typename T, uint N>
struct Cells
{
    std::vector<std::array<uint, static_cast<uint>(1) << N>> edgeIndices;
    std::vector<Vec<T, N>> cellAverages;
    std::vector<T> cellVolumes;
    std::vector<Vec<T, N>> cellCentersPositions;
};
