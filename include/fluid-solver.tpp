#pragma once
#include <algorithm>
#include <iostream>
// #include <SFML/Window/Keyboard.hpp>

// #include "fluid-solver.hpp"

template <typename T, uint N>
Solver<T, N>::Solver(Edges<T, N>& _edges, Cells<T, N>& _cells, T _a)
    : edges(_edges), cells(_cells)
{
    dimension = N;
    numberOfCells = cells.cellAverages.size();
    numberOfEdges = edges.edgeProperty.size();
    assert(cells.cellVolumes.size() == numberOfCells);
    if(verifyDomainIntegrity())
    {
        std::cout << "Domain integrity verified." << std::endl;
    }
    else
    {
        std::cerr << "Domain integrity verification failed." << std::endl;
        exit(EXIT_FAILURE);
    }
    a = _a;
}


template <typename T, uint N>
bool Solver<T,N>::verifyDomainIntegrity() const
{
    std::cout << "Verifying domain integrity..." << std::endl;
    constexpr uint n1 = 1 << N;
    constexpr uint n2 = 1 << (N - 1);
    std::cout<< "Number of cells: " << numberOfCells << std::endl;
    std::cout<< "Number of edges: " << numberOfEdges << std::endl;
    std::cout<< "Number of edge indices "<< cells.edgeIndices.size() << std::endl;
    for (uint i = 1; i < numberOfEdges-1; i++)
    {
        // std::vector<std::array<Vec<T, N>, static_cast<uint>(1) << (N - 1)>>
        EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[i];
        std::array<Vec<T, N>, n2> vertices = edges.vertices[i];
        
        uint cell1Index = edgeProperty.connectingCells.first;
        uint cell2Index = edgeProperty.connectingCells.second;
        if (cell1Index >= cells.edgeIndices.size())
        {
            std::cerr << "Invalid cell1Index: " << cell1Index
                    << ", number of cells: " << cells.edgeIndices.size()
                    << '\n';
            return false;
        }

        if (cell2Index >= cells.edgeIndices.size())
        {
            std::cerr << "Invalid cell2Index: " << cell2Index
                    << ", number of cells: " << cells.edgeIndices.size()
                    << '\n';
            return false;
        }
        std::array<uint, n1>& cell1EdgeIndices = cells.edgeIndices[cell1Index];
        Vec<T, N> cell1Center = cells.cellCentersPositions[cell1Index];
        std::array<uint, n1>& cell2EdgeIndices = cells.edgeIndices[cell2Index];
        Vec<T, N> cell2Center = cells.cellCentersPositions[cell2Index];
        std::cout<<"Beginning checks for edge "<<i<<" connecting cells "<<cell1Index<<" and "<<cell2Index<<std::endl;
        if (std::find(cell1EdgeIndices.begin(), cell1EdgeIndices.end(), i) == cell1EdgeIndices.end())
        {
            std::cerr << "Error: Edge " << i << " is not found in the edge indices of cell " << cell1Index << std::endl;
            std::cerr << "Edge indices of cell " << cell1Index << ": ";
            for (const auto& edgeIndex : cell1EdgeIndices)
            {
                std::cerr << edgeIndex << " ";
            }
            std::cerr << std::endl;
            return false;
        }
        if (std::find(cell2EdgeIndices.begin(), cell2EdgeIndices.end(), i) == cell2EdgeIndices.end())
        {
            std::cerr << "Error: Edge " << i << " is not found in the edge indices of cell " << cell2Index << std::endl;
            std::cerr << "Edge indices of cell " << cell2Index << ": ";
            for (const auto& edgeIndex : cell2EdgeIndices)
            {
                std::cerr << edgeIndex << " ";
            }
            std::cerr << std::endl;
            return false;
        }
        std::cout<<"Edge "<<i<<" is found in the edge indices of both cells."<<std::endl;
        // Okay check one done. now for checking if the edge normal is pointing from cell1 to cell2.
        for(uint vertexId = 0; vertexId < vertices.size(); vertexId++)
        {
            Vec<T, N> vertex = vertices[vertexId];
            Vec<T, N> cell1ToVertex = vertex - cell1Center;
            Vec<T, N> cell2ToVertex = vertex - cell2Center;
            T dotProduct1 = 0;
            T dotProduct2 = 0;
            if(i== numberOfCells-1) continue;
            for(uint dim = 0; dim < N; dim++)
            {
                dotProduct1 += edgeProperty.normal[dim] * cell1ToVertex[dim];
                dotProduct2 += edgeProperty.normal[dim] * cell2ToVertex[dim];
            }
            if(dotProduct1 <= 0)
            {
                std::cerr << "Error: Edge " << i << " normal is pointing towards cell " << cell1Index << std::endl;
                return false;
            }
            if(dotProduct2 >= 0)
            {
                std::cerr << "Error: Edge " << i << " normal is away from cell " << cell2Index << std::endl;
                return false;
            }
        }
    }
    return true;
}


template <typename T, uint N>
void Solver<T, N>::calculateFluxes()
{
    for (uint i = 0; i < numberOfEdges; i++)
    {
        EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[i];
        Vec<T, N>& cell1Average = cells.cellAverages[edgeProperty.connectingCells.first];
        Vec<T, N>& cell2Average = cells.cellAverages[edgeProperty.connectingCells.second];

        T velFlux1 = 0;
        T velFlux2 = 0;

        // Calculate fluxes based on the edge normal and cell averages

        for (uint j = 0; j < N; j++)
        {
            velFlux1 += edgeProperty.normal[j] * cell1Average[j];
            velFlux2 += edgeProperty.normal[j] * cell2Average[j]; // First order for now
        }

        edgeProperty.normalVelFlux[0] = velFlux1;
        edgeProperty.normalVelFlux[1] = velFlux2;
        // edgeProperty.numericalFlux = (velFlux1 * (a + absVal(a)))* 0.5 + (velFlux2 * (a - absVal(a)))* 0.5; // Rusanov flux for now
        edgeProperty.numericalFlux = (a * (edgeProperty.normalVelFlux[0] + edgeProperty.normalVelFlux[1])  - absVal(a)*(edgeProperty.normalVelFlux[1] - edgeProperty.normalVelFlux[0])) * 0.5;

    }

}
template <typename T, uint N>
void Solver<T, N>::applyBoundaryConditions()
{
    // Apply periodic boundary conditions for now
    // lastEdge.normalVel
    // EdgeProperty<T, N>& firstEdge = edges.edgeProperty[0];
    // EdgeProperty<T, N>& lastEdge = edges.edgeProperty[numberOfEdges - 1];
    // T averageFlux = (firstEdge.numericalFlux + lastEdge.numericalFlux) * 0.5;
    // firstEdge.numericalFlux = averageFlux;
    // lastEdge.numericalFlux = averageFlux;
    // lastEdge.numericalFlux = firstEdge.numericalFlux;

    cells.cellAverages[numberOfCells-1][0] = cells.cellAverages[numberOfCells - 2][0]; 
    edges.edgeProperty[numberOfEdges-1].numericalFlux = edges.edgeProperty[numberOfEdges - 2].numericalFlux;
    cells.cellAverages[0][0] = cells.cellAverages[1][0];
    edges.edgeProperty[0].numericalFlux = edges.edgeProperty[1].numericalFlux;
}

template <typename T, uint N>
void Solver<T, N>::updateCellAverages(T dt){
    for(uint i = 0; i < numberOfCells; i++)
    {
        std::array<uint, static_cast<uint>((1u)<<N)>& cellEdgeIndices = cells.edgeIndices[i];
        // T& cellAverage = cells.cellAverages[i];
        const T cellVolume = cells.cellVolumes[i];

        T fluxLeaving = 0;
        for (uint j = 0; j < cellEdgeIndices.size(); j++)
        {
            uint edgeIndex = cellEdgeIndices[j];
            EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[edgeIndex];
            fluxLeaving += edgeProperty.numericalFlux * (static_cast<T>((edgeProperty.connectingCells.first == i))*2 - static_cast<T>(1)); // Branchless version of the above if-else statement. If fluxExiting is true, multiply by 1, else multiply by -1.
        }
        cells.cellAverages[i][0] -= fluxLeaving * dt / cellVolume;
    }
}


