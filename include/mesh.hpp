#pragma once

#include <memory>

#include "constants.hpp"
#include "vec.hpp"
#include "cells.hpp"
#include "domains.hpp"


enum class MeshType
{
    NONE, 
    CARTESIAN,
    UNSTRUCTURED,
};

template <typename T, uint N>
struct Mesh{
    uint numCells = 0;
    uint numEdges = 0;
    MeshType meshType = MeshType::CARTESIAN;
    Edges<T, N> edges;
    std::pair<Vec<T, N>, Vec<T, N>> boundingBox;
    // std::vector<Index> boundaryEdgesIndices;
    // std::vector<Index> boundaryCellsIndices;
    // std::vector<Index> domainEdgesIndices;
    Cells<T, N> cells;
    void fixEdgesConnectivityOrder();
    void identifyBoundaryEdgesAndCells(const std::vector<std::pair<uint, std::unique_ptr<Domain<T, N>>>>& boundaryDomains);
};

class MeshGenerator{
    uint dimensions;
    uint numberOfCells;
    MeshType meshType;
    Vec<double, 2> domainBounds; // xMin, xMax
    public:
    MeshGenerator(uint _dimensions, Vec<double, 2> _domainBounds, uint _numberOfCells, MeshType _meshType = MeshType::CARTESIAN);
    Mesh<double, 1> generateMesh();
};


