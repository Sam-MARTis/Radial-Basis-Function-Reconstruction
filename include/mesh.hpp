#pragma once

#include "constants.hpp"
#include "vec.hpp"
#include "cells.hpp"

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
    Cells<T, N> cells;
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


