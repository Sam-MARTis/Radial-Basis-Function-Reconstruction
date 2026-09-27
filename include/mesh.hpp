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
    uint numBoundaries = 1;
    Edges<T, N> edges;
    std::pair<Vec<T, N>, Vec<T, N>> boundingBox;
    std::vector<Index> boundaryEdgesIndices;
    std::vector<Index> boundaryCellsIndices;
    std::vector<Index> domainEdgesIndices;
    Cells<T, N> cells;
    std::vector<std::pair<uint, std::vector<Index>>> boundaryIdentifiers;
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
        // neighbourhood.neighbourCellIndices.clear();
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
        neighbourhood.rbfCoefficients = std::array<std::vector<T>, N+2>();
        mesh.cells.neighbourhoods.push_back(neighbourhood);
    }
}


template<typename T, uint N>
void Mesh<T, N>::fixEdgesConnectivityOrder()
{
    // std::cout << "fixEdgesConnectivityOrder not implemented yet." << std::endl;

    for (Index edgeIdx = 0; edgeIdx < numEdges; edgeIdx++)
    {
        EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[edgeIdx];
        const Index sourceCellIndex = edgeProperty.connectingCells.first;
        const Index sinkCellIndex = edgeProperty.connectingCells.second;
        std::array<T, N>& edgeNormal = edgeProperty.normal;
        if (sourceCellIndex >= numCells || sinkCellIndex >= numCells)
        {
            assert(((sourceCellIndex>= numCells) && (sinkCellIndex >= numCells))==0);
            Vec<T, N> edgeCenter(0);
            for (uint i=0; i<(1<<(N-1)); i++)
            {
                edgeCenter += edges.vertices[edgeIdx][i];
            }
            // auto edgeVertex1 = edges.vertices[edgeIdx][0];
            edgeProperty.measure = (edges.vertices[edgeIdx][1] - edges.vertices[edgeIdx][0]).norm();
            // std::cout<<"Edge "<<edgeIdx<<" measure: "<<edgeProperty.measure<<std::endl;
            edgeCenter *= static_cast<T>(1.0/(1<<(N-1)));
            if (sourceCellIndex>= numCells)
            {
                const Vec<T, N>& sinkCellCenter = cells.cellCentersPositions[sinkCellIndex];
                const Vec<T, N> directionVec = sinkCellCenter - edgeCenter;
                T dotVal = Vec<T, N>::dot(directionVec, edgeNormal);
                if (dotVal < 0) edgeProperty.connectingCells = {sinkCellIndex, sourceCellIndex};
            }
            if (sinkCellIndex >= numCells)
            {
                const Vec<T, N>& sourceCellCenter = cells.cellCentersPositions[sourceCellIndex];
                const Vec<T, N> directionVec = edgeCenter - sourceCellCenter;
                T dotVal = Vec<T, N>::dot(directionVec, edgeNormal);
                if (dotVal < 0) edgeProperty.connectingCells = {sinkCellIndex, sourceCellIndex};

            }

            /*
             * This branch implies the edge is a booundary edge with only one connecting cell.
             * In the future we might wanna do some aligning but for now that is not necessary.
             * Fluxes for boundary and stuff will be handled separately in the constraints
             *
             * Note from future Sep8 2026, the above was implemented and normals were aligned
             * source sink relation with normal is followed for boundary edges as well
             */
            continue;
        }
        const Vec<T, N>& sourceCellCenter = cells.cellCentersPositions[sourceCellIndex];
        const Vec<T, N>& sinkCellCenter = cells.cellCentersPositions[sinkCellIndex];
        const Vec<T, N> directionVec = sinkCellCenter - sourceCellCenter;
        T dotVal = Vec<T, N>::dot(directionVec, edgeNormal);
        // for (Index i = 0; i < N; i++) dotVal += directionVec[i] * edgeNormal[i];
        if (dotVal < 0)
        {
            edgeProperty.connectingCells = {sinkCellIndex, sourceCellIndex};
            // for (uint edgeNormalIdx = 0; edgeNormalIdx < N; edgeNormalIdx++)
            // {
            //     edgeNormal[edgeNormalIdx] *= -1.0;
            // }
        }
    }
}




template<typename T, uint N>
void Mesh<T, N>::identifyBoundaryEdgesAndCells(const  std::vector<std::pair<uint, std::unique_ptr<Domain<T, N>>>>& boundaryDomains)
{
    const uint numSuppliedBoundaryDomains = boundaryDomains.size();
    numBoundaries = numSuppliedBoundaryDomains;
    std::unordered_map<uint, uint> boundaryIdToIndexMapping;
    for (Index i=0; i<numSuppliedBoundaryDomains; i++)
    {

        const auto& [boundaryId, boundaryDomain] = boundaryDomains[i];
        boundaryIdentifiers.push_back({boundaryId, {}});
        boundaryIdToIndexMapping.insert({boundaryId, i});
    }
    // const Index numEdges = this->numEdges;
    for (Index edgeIdx = 0; edgeIdx < numEdges; edgeIdx++)
    {
        const EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[edgeIdx];
        const Index sourceCellIndex = edgeProperty.connectingCells.first;
        const Index sinkCellIndex = edgeProperty.connectingCells.second;
        if (sourceCellIndex >= numCells || sinkCellIndex >= numCells)
        {
            if (sourceCellIndex>=numCells)
            {
                boundaryCellsIndices.push_back(sinkCellIndex);
            }
            else
            {
                boundaryCellsIndices.push_back(sourceCellIndex);
            }
            boundaryEdgesIndices.push_back(edgeIdx);
            for (Index boundaryCheckIdx = 0; boundaryCheckIdx < numSuppliedBoundaryDomains; boundaryCheckIdx++)
            {
                const auto& [boundaryIdTag, boundaryDomain] = boundaryDomains[boundaryCheckIdx];
                const Vec<T, N>& edgeVertex1 = edges.vertices[edgeIdx][0];
                const Vec<T, N>& edgeVertex2 = edges.vertices[edgeIdx][1];
                if (boundaryDomain->isInside(edgeVertex1) && boundaryDomain->isInside(edgeVertex2))
                {
                    uint boundaryIdx = boundaryIdToIndexMapping.at(boundaryIdTag);
                    boundaryIdentifiers[boundaryIdx].second.push_back(edgeIdx);
                    break;
                }
            }
        }else
        {
            domainEdgesIndices.push_back(edgeIdx);
        }
    }

}



