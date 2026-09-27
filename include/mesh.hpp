#pragma once

#include <memory>
#include <set>

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
    static void initializeNeighbourhood(Mesh<T, N>& mesh, uint depth);
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


template <typename T, uint N>
void Mesh<T, N>::initializeNeighbourhood(Mesh<T, N>& mesh, uint depth)
{
    assert(mesh.cells.neighbourhoods.size() == 0);
    for (Index cellId = 0; cellId < mesh.numCells; cellId++)
    {
        Neighbourhood<T, N> neighbourhood;
        neighbourhood.numNeighbours = 0;
        neighbourhood.neighbourCellIndices.clear();
        neighbourhood.rbfCoefficients.clear();
        neighbourhood.ALocalInv = Matrix<T>(depth, depth);
        std::set<Index> neighbourhoodSets = {cellId};
        for (uint currentDepth = 0; currentDepth < depth; currentDepth++)
        {
            std::vector<Index> thisDepthNeighbours;
            for (const Index& neighbourId : neighbourhoodSets)
            {
                std::vector<Index> cellConnectingEdges = mesh.cells.edgeIndices[neighbourId];
                for (const Index& edgeId : cellConnectingEdges)
                {
                    EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[edgeId];
                    Index neighbourCellId = (edgeProperty.connectingCells.first == neighbourId) ? edgeProperty.connectingCells.second : edgeProperty.connectingCells.first;
                    if (neighbourCellId < mesh.numCells)
                    {
                        thisDepthNeighbours.push_back(neighbourCellId);
                        // thisDepthNeighbours.insert(neighbourCellId);
                    }
                }
            }
            for (const Index& neighbourId : thisDepthNeighbours)
            {
                neighbourhoodSets.insert(neighbourId);
            }
        }
        // neighbourhoodSets.erase(cellId);
        for (const Index& neighbourId : neighbourhoodSets)
        {
            neighbourhood.neighbourCellIndices.push_back(neighbourId);
        }
        neighbourhood.numNeighbours = neighbourhood.neighbourCellIndices.size();
    }
}

