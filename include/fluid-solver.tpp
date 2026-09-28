#pragma once
#include <algorithm>
#include <iostream>
#pragma <omp.h>
// #include <SFML/Window/Keyboard.hpp>

// #include "fluid-solver.hpp"

// template <typename T, uint N>
// Solver<T, N>::Solver(Edges<T, N>& _edges, Cells<T, N>& _cells, T _a)
//     : edges(_edges), cells(_cells)
// {
//     dimension = N;
//     numberOfCells = cells.cellAverages.size();
//     numberOfEdges = edges.edgeProperty.size();
//     assert(cells.cellVolumes.size() == numberOfCells);
//     // if(verifyDomainIntegrity())
//     // {
//     //     std::cout << "Domain integrity verified." << std::endl;
//     // }
//     // else
//     // {
//     //     std::cerr << "Domain integrity verification failed." << std::endl;
//     //     exit(EXIT_FAILURE);
//     // }
//     a = _a;
// }

template <typename T, uint N>
Solver<T, N>::Solver(Mesh<T, N>& mesh)
    : mesh(mesh)
{
    std::cout << "Initializing Solver..." << std::endl;
    dimension = N;
    numberOfCells = mesh.cells.cellAverages.size();
    numberOfEdges = mesh.edges.edgeProperty.size();
    assert(mesh.cells.cellVolumes.size() == numberOfCells);
    boundaryConditions.resize(mesh.numBoundaries);
}

//
// template <typename T, uint N>
// bool Solver<T,N>::verifyDomainIntegrity() const
// {
//     std::cout << "Verifying domain integrity..." << std::endl;
//     constexpr uint n1 = 1 << N;
//     constexpr uint n2 = 1 << (N - 1);
//     std::cout<< "Number of cells: " << numberOfCells << std::endl;
//     std::cout<< "Number of edges: " << numberOfEdges << std::endl;
//     std::cout<< "Number of edge indices "<< cells.edgeIndices.size() << std::endl;
//     for (uint i = 1; i < numberOfEdges-1; i++)
//     {
//         // std::vector<std::array<Vec<T, N>, static_cast<uint>(1) << (N - 1)>>
//         EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[i];
//         std::array<Vec<T, N>, n2> vertices = edges.vertices[i];
//
//         uint cell1Index = edgeProperty.connectingCells.first;
//         uint cell2Index = edgeProperty.connectingCells.second;
//         if (cell1Index >= cells.edgeIndices.size())
//         {
//             std::cerr << "Invalid cell1Index: " << cell1Index
//                     << ", number of cells: " << cells.edgeIndices.size()
//                     << '\n';
//             return false;
//         }
//
//         if (cell2Index >= cells.edgeIndices.size())
//         {
//             std::cerr << "Invalid cell2Index: " << cell2Index
//                     << ", number of cells: " << cells.edgeIndices.size()
//                     << '\n';
//             return false;
//         }
//         std::array<uint, n1>& cell1EdgeIndices = cells.edgeIndices[cell1Index];
//         Vec<T, N> cell1Center = cells.cellCentersPositions[cell1Index];
//         std::array<uint, n1>& cell2EdgeIndices = cells.edgeIndices[cell2Index];
//         Vec<T, N> cell2Center = cells.cellCentersPositions[cell2Index];
//         std::cout<<"Beginning checks for edge "<<i<<" connecting cells "<<cell1Index<<" and "<<cell2Index<<std::endl;
//         if (std::find(cell1EdgeIndices.begin(), cell1EdgeIndices.end(), i) == cell1EdgeIndices.end())
//         {
//             std::cerr << "Error: Edge " << i << " is not found in the edge indices of cell " << cell1Index << std::endl;
//             std::cerr << "Edge indices of cell " << cell1Index << ": ";
//             for (const auto& edgeIndex : cell1EdgeIndices)
//             {
//                 std::cerr << edgeIndex << " ";
//             }
//             std::cerr << std::endl;
//             return false;
//         }
//         if (std::find(cell2EdgeIndices.begin(), cell2EdgeIndices.end(), i) == cell2EdgeIndices.end())
//         {
//             std::cerr << "Error: Edge " << i << " is not found in the edge indices of cell " << cell2Index << std::endl;
//             std::cerr << "Edge indices of cell " << cell2Index << ": ";
//             for (const auto& edgeIndex : cell2EdgeIndices)
//             {
//                 std::cerr << edgeIndex << " ";
//             }
//             std::cerr << std::endl;
//             return false;
//         }
//         std::cout<<"Edge "<<i<<" is found in the edge indices of both cells."<<std::endl;
//         // Okay check one done. now for checking if the edge normal is pointing from cell1 to cell2.
//         for(uint vertexId = 0; vertexId < vertices.size(); vertexId++)
//         {
//             Vec<T, N> vertex = vertices[vertexId];
//             Vec<T, N> cell1ToVertex = vertex - cell1Center;
//             Vec<T, N> cell2ToVertex = vertex - cell2Center;
//             T dotProduct1 = 0;
//             T dotProduct2 = 0;
//             if(i== numberOfCells-1) continue;
//             for(uint dim = 0; dim < N; dim++)
//             {
//                 dotProduct1 += edgeProperty.normal[dim] * cell1ToVertex[dim];
//                 dotProduct2 += edgeProperty.normal[dim] * cell2ToVertex[dim];
//             }
//             if(dotProduct1 <= 0)
//             {
//                 std::cerr << "Error: Edge " << i << " normal is pointing towards cell " << cell1Index << std::endl;
//                 return false;
//             }
//             if(dotProduct2 >= 0)
//             {
//                 std::cerr << "Error: Edge " << i << " normal is away from cell " << cell2Index << std::endl;
//                 return false;
//             }
//         }
//     }
//     return true;
// }



template <typename T, uint N>
void Solver<T, N>::calculateFluxes()
{
    // assert(N==2);
    #pragma omp parallel for
    for (const uint domainEdgeIndex: mesh.domainEdgesIndices)
    {
        EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[domainEdgeIndex];        // if (edgeProperty.connectingCells.second >= mesh.numCells)
        const Vec<T, N+2> UL = edgeProperty.reconstructedPushedU.first;
        const Vec<T, N+2> UR = edgeProperty.reconstructedPushedU.second;
        const std::array<T, N>& edgeNormal = edgeProperty.normal;
        std::array<std::array<T, N>, N> directionVectors;
        directionVectors[0] = edgeNormal;
        if (N==2)
        {
            directionVectors[1] = {-edgeNormal[1], edgeNormal[0]};
        }
        else
        {
            std::cerr << "Error: Unsupported dimension N=" << N << " for flux computation." << std::endl;
            exit(EXIT_FAILURE);
        }
        {
            // Profiler profiler("calculateFluxes::RoeFlux");
            edgeProperty.numericalFlux = RoeFlux<T, N>::computeNumericalFlux(UL, UR, directionVectors);
        }
            // edgeProperty.numericalFlux = RussanovFlux<T, N>::computeNumericalFlux(cellLAverage, cellRAverage, directionVectors);
        /*
         * Put the task of computing out vs in on the cell update function.
         * Remember, normal points from connectingCells.first to connectingCells.second.
         */
    }

}


// template <typename T, uint N>
// void Solver<T, N>::calculateRBFBasedFluxes(const RBFs<T, N>& rbfs, const T dt)
// {
//     const bool aRight = a > 0;
//     for (uint i=0; i<numberOfEdges; i++)
//     {
//         EdgeProperty<T, N>& edgeProperty = edges.edgeProperty[i];
//         const T edgePosition = edges.vertices[i][0][0];
//         const T lb = a>0? edgePosition - dt * a : edgePosition;
//         const T ub = a<0? edgePosition - dt * a : edgePosition;
//         // const T lb = edgePosition - dt*a*(static_cast<T>(aRight));
//         // const T ub = edgePosition - dt*a*(static_cast<T>(!aRight));
//         assert(lb<ub);
//         // std::cout<< "ub: " << ub << " lb: " << lb << std::endl;
//         // std::cout << "absVal(ub - lb)" << absVal(ub - lb) << std::endl;
//         assert(absVal(ub - lb) > EPSILON1);
//         const T rbfValue = rbfs.integrate(lb, ub);
//         edgeProperty.numericalFlux = rbfValue;
//         // std::cout <<edgeProperty.numericalFlux;
//     }
// }


template <typename T, uint N>
void Solver<T, N>::applyBoundaryConditions()
{
    for (uint boundaryIdx=0; boundaryIdx<mesh.numBoundaries; boundaryIdx++)
    {
        switch (boundaryConditions[boundaryIdx].first)
        {
            case BOUNDARY_CONSTRAINT::DIRICHLET:
                Constraint<T, N>::DirichletConditionAtEdges(mesh, mesh.boundaryIdentifiers[boundaryIdx].second, boundaryConditions[boundaryIdx].second);
                break;
            case BOUNDARY_CONSTRAINT::EXTRAPOLATED:
                Constraint<T, N>::ExtrapolatedFluxAtEdges(mesh, mesh.boundaryIdentifiers[boundaryIdx].second);
                break;
            case BOUNDARY_CONSTRAINT::SOLID_WALL:
                Constraint<T, N>::SolidWallConditionAtEdges(mesh, mesh.boundaryIdentifiers[boundaryIdx].second);
                break;
            default:
                std::cerr << "Error: Unsupported boundary condition type." << std::endl;
                exit(EXIT_FAILURE);
        }
    }
}

template <typename T, uint N>
void Solver<T, N>::updateCellAverages(T dt){
#pragma omp parallel for
    for(Index cellIdx = 0; cellIdx < numberOfCells; cellIdx++)
    {
        std::vector<Index>& cellEdgeIndices = mesh.cells.edgeIndices[cellIdx];
        // T& cellAverage = mesh.cells.cellAverages[i];
        const T invCellVolume = 1.0/mesh.cells.cellVolumes[cellIdx];

        Vec<T, N+2> fluxLeaving(0);
        // std::cout<<"Edge 0 measure: "<<mesh.edges.edgeProperty[cellEdgeIndices[0]].measure<<std::endl;
        for (unsigned int edgeIndex : cellEdgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[edgeIndex];
            // std::cout<<"Measure of edge is: "<<edgeProperty.measure<<std::endl;
            fluxLeaving += edgeProperty.numericalFlux * edgeProperty.measure * (static_cast<T>((edgeProperty.connectingCells.first == cellIdx)) * 2 - static_cast<T>(1)); // Branchless version of the above if-else statement. If fluxExiting is true, multiply by 1, else multiply by -1.
        }
        const T rhoBeforeUpdate = mesh.cells.cellAverages[cellIdx][0];
        assert(rhoBeforeUpdate > 0);
        // assert(fluxLeaving)

        mesh.cells.cellAverages[cellIdx] -= fluxLeaving * dt*invCellVolume;
        if (!(mesh.cells.cellAverages[cellIdx][0]>0))
        {
            const auto U = mesh.cells.cellAverages[cellIdx];
            std::cout<<"ERROR: Rho: "<< U[0] <<", E: "<< U[N-1] <<std::endl;
        }
        assert(mesh.cells.cellAverages[cellIdx][0] > 0);
        for (uint fieldIdx = 0; fieldIdx<N+2; fieldIdx++)
        {
            assert(!std::isnan(mesh.cells.cellAverages[cellIdx][fieldIdx]));
        }
    }
}

template <typename T, uint N>
void Solver<T, N>::updateNeighbourhoodRBFCoeffsAndExtremums()
{
#pragma omp parallel for
    for (Index cellId = 0; cellId<mesh.numCells; cellId++)
    {
        Neighbourhood<T, N>& neighbourhood = mesh.cells.neighbourhoods[cellId];
        const std::vector<Index>& neighbourCellIndices = neighbourhood.neighbourCellIndices;
        const uint numNeighbours = neighbourhood.numNeighbours;
        std::vector<Vec<T, N+2>> localNeighbourAverages(numNeighbours);
        Vec<T, N+2>& minU = neighbourhood.extremeValues.first;
        Vec<T, N+2>& maxU = neighbourhood.extremeValues.second;
        for (uint fieldIdx = 0; fieldIdx < N+2; fieldIdx++) {
            minU[fieldIdx] = std::numeric_limits<T>::max();
            maxU[fieldIdx] = std::numeric_limits<T>::lowest();
        }
        for (uint neighbourCellLocalIdx=0; neighbourCellLocalIdx<numNeighbours; neighbourCellLocalIdx++)
        {
            const Index neighbourCellGlobalIdx = neighbourCellIndices[neighbourCellLocalIdx];
            const Vec<T,N+2>& UNeighbour = mesh.cells.cellAverages[neighbourCellGlobalIdx];
            localNeighbourAverages[neighbourCellLocalIdx] = mesh.cells.cellAverages[neighbourCellGlobalIdx];
            for (uint fieldIdx=0; fieldIdx<N+2; fieldIdx++)
            {
                if (UNeighbour[fieldIdx] < minU[fieldIdx]) minU[fieldIdx] = UNeighbour[fieldIdx];
                if (UNeighbour[fieldIdx] > maxU[fieldIdx]) maxU[fieldIdx] = UNeighbour[fieldIdx];
            }
        }

        // std::array<std::vector<T>, N+2> rbfCoefficients;
        for (uint fieldVarIdx=0; fieldVarIdx<N+2; fieldVarIdx++)
        {
            std::vector<T> coeffs(numNeighbours);
            for (uint neighbourCellLocalIdx=0; neighbourCellLocalIdx<numNeighbours; neighbourCellLocalIdx++)
            {
                coeffs[neighbourCellLocalIdx] = localNeighbourAverages[neighbourCellLocalIdx][fieldVarIdx];
            }
            neighbourhood.rbfCoefficients[fieldVarIdx] = MatrixSolver<T>::matrixStdVectorMultiply(neighbourhood.ALocalInv, coeffs);
        }
        // mesh.cells.neighbourhoods[cellId].rbfCoefficients = rbfCoefficients;
        // mesh.cells.neighbourhoods[cellId].rbfCoefficients = rbfCoefficients;
    }
}
/*

template <typename T, uint N>
void Solver<T, N>::reconstructSolutionAtEdges()
{
    #pragma omp parallel for
    for (Index cellIdx = 0; cellIdx < numberOfCells; cellIdx++)
    {
        const Neighbourhood<T, N>& neighbourhood = mesh.cells.neighbourhoods[cellIdx];
        const Vec<T, N+2>& Uj = mesh.cells.cellAverages[cellIdx];
        const std::vector<Index>& cellEdgeIndices = mesh.cells.edgeIndices[cellIdx];
        // const uint numNeighbours = neighbourhood.neighbourCellsIndices.size();
        const Vec<T, N>& cellCenter = mesh.cells.cellCentersPositions[cellIdx];
        const Vec<Vec<T, N>, N+2>& cellGradients = mesh.cells.reconstructedCellGradients[cellIdx];
        for (Index edgeIdx : cellEdgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[edgeIdx];
            const Vec<T, N>& edgeEvaluationPoint = mesh.edges.edgeEvaluationPoints[edgeIdx];
            const Vec<T, N> dr = edgeEvaluationPoint - cellCenter;
            const Vec<T, N+2> dU = VecMath::matrixMultiply<T, N+2, N>(cellGradients, dr);
            const Vec<T, N+2> Ui = Uj + dU;
            const T epsilonSq = LIMITER_K * LIMITER_K * edgeProperty.measure * edgeProperty.measure;
            const Vec<T, N+2> psi = Limiters<T, N+2>::Venkatakrishnan(Uj, Ui, neighbourhood.extremeValues.first, neighbourhood.extremeValues.second, epsilonSq);
            const Vec<T, N+2> Ulimited = Uj + VecMath::componentwiseMultiply<T, N+2>(psi, dU);
            if (edgeProperty.connectingCells.first == cellIdx) {
                edgeProperty.reconstructedPushedU.first = Ulimited;
            } else if (edgeProperty.connectingCells.second == cellIdx) {
                edgeProperty.reconstructedPushedU.second = Ulimited;
            } else {
                assert(false && "edge does not belong to this cell");
            }
        }
    }
}
 */


template <typename T, uint N>
void Solver<T, N>::prescribeBoundaryConditions(uint boundaryIndex, BOUNDARY_CONSTRAINT boundaryType)
{
    assert(boundaryType != BOUNDARY_CONSTRAINT::DIRICHLET);
    prescribeBoundaryConditions(boundaryIndex, boundaryType, Vec<T, N+2>(0));
}

template <typename T, uint N>
void Solver<T, N>::prescribeBoundaryConditions(uint boundaryIndex, BOUNDARY_CONSTRAINT boundaryType, Vec<T, N+2> prescribedBoundaryValue)
{
    // assert(boundaryType != BOUNDARY_CONSTRAINT::DIRICHLET);
    assert(boundaryConditions.size() == mesh.numBoundaries);
    for (uint i=0; i<mesh.numBoundaries; i++)
    {
        if (mesh.boundaryIdentifiers[i].first == boundaryIndex)
        {
            boundaryConditions[i] = {boundaryType, prescribedBoundaryValue};
        }
    }
}


template <typename T, uint N>
void Solver<T, N>::reconstructSolutionAtEdgesViaRBFs(const RBFs<T, N>& rbfs)
{
    for (Index currentCellIdx = 0; currentCellIdx < numberOfCells; currentCellIdx++)
    {
        const Neighbourhood<T, N>& neighbourhood = mesh.cells.neighbourhoods[currentCellIdx];
        const Vec<T, N+2>& Uj = mesh.cells.cellAverages[currentCellIdx];

        const std::vector<Index>& cellEdgeIndices = mesh.cells.edgeIndices[currentCellIdx];
        const Vec<T, N>& cellCenter = mesh.cells.cellCentersPositions[currentCellIdx];
        // const Vec<Vec<T, N>, N+2>& cellGradients = mesh.cells.reconstructedCellGradients[currentCellIdx];
        for (Index edgeIdx : cellEdgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[edgeIdx];
            const Vec<T, N>& edgeEvaluationPoint = mesh.edges.edgeEvaluationPoints[edgeIdx];
            Vec<T, N+2> Ui(0);
            for (Index fieldIdx = 0; fieldIdx < N+2; fieldIdx++)
            {
                Ui[fieldIdx] = rbfs.valueLocal(edgeEvaluationPoint, neighbourhood.neighbourCellIndices, neighbourhood.rbfCoefficients[fieldIdx]);
            }
            const Vec<T, N+2> dU = Ui - Uj;
            const T epsilonSq = LIMITER_K * LIMITER_K * edgeProperty.measure * edgeProperty.measure;
            const Vec<T, N+2> psi = Limiters<T, N+2>::Venkatakrishnan(Uj, Ui, neighbourhood.extremeValues.first, neighbourhood.extremeValues.second, epsilonSq);
            const Vec<T, N+2> Ulimited = Uj + VecMath::componentwiseMultiply<T, N+2>(psi, dU);
            if (edgeProperty.connectingCells.first == currentCellIdx) {
                edgeProperty.reconstructedPushedU.first = Ulimited;
            } else if (edgeProperty.connectingCells.second == currentCellIdx) {
                edgeProperty.reconstructedPushedU.second = Ulimited;
            } else {
                assert(false && "edge does not belong to this cell");
            }
        }
    }
}

template <typename T, uint N>
void Solver<T, N>::reconstructSolutionAtEdgesPiecewise()
{
    for (Index currentCellIdx = 0; currentCellIdx < numberOfCells; currentCellIdx++)
    {
        const Vec<T, N+2>& Uj = mesh.cells.cellAverages[currentCellIdx];
        const std::vector<Index>& cellEdgeIndices = mesh.cells.edgeIndices[currentCellIdx];
        for (Index edgeIdx : cellEdgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[edgeIdx];
            if (edgeProperty.connectingCells.first == currentCellIdx) {
                edgeProperty.reconstructedPushedU.first = Uj;
            } else if (edgeProperty.connectingCells.second == currentCellIdx) {
                edgeProperty.reconstructedPushedU.second = Uj;
            } else {
                assert(false && "edge does not belong to this cell");
            }
        }
    }
}
