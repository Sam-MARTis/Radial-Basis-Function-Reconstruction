#pragma once
#include "vec.hpp"
// #include "mesh.hpp"
#include "domains.hpp"
#include "physics.hpp"


enum class BOUNDARY_CONSTRAINT
{
    EXTRAPOLATED,
    DIRICHLET,
    SOLID_WALL,
};

template <typename T, uint N>
class Constraint
{
public:
    static void DirichletCondition(Mesh<T, N>& mesh, const Domain<T, N>& region, const Vec<T, N+2>& value)
    {
        for (Index i=0; i<mesh.numCells; i++)
        {
            const Vec<T, N>& cellCenter = mesh.cells.cellCentersPositions[i];
            if (region.isInside(cellCenter))
            {
                mesh.cells.cellAverages[i] = value;
            }
        }
    }
    static void DirichletConditionAtEdges(Mesh<T, N>& mesh, const std::vector<Index>& edgeIndices, const Vec<T, N+2>& UPrescribed)
    {
        for (const Index boundaryEdgeIndex : edgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[boundaryEdgeIndex];
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
            if (edgeProperty.connectingCells.first >= mesh.numCells)
            {
                // const Vec<T, N+2>& UR = mesh.cells.cellAverages[edgeProperty.connectingCells.second];
                const Vec<T, N+2>& UR = edgeProperty.reconstructedPushedU.second;
                edgeProperty.numericalFlux = RoeFlux<T, N>::computeNumericalFlux(UPrescribed, UR, directionVectors);
            }
            else if (edgeProperty.connectingCells.second >= mesh.numCells)
            {
                // const Vec<T, N+2>& UL = mesh.cells.cellAverages[edgeProperty.connectingCells.first];
                const Vec<T, N+2>& UL = edgeProperty.reconstructedPushedU.first;
                edgeProperty.numericalFlux = RoeFlux<T, N>::computeNumericalFlux(UL, UPrescribed, directionVectors);
            }

        }
    }
    static void ExtrapolatedFluxAtEdges(Mesh<T, N>& mesh, const std::vector<Index>& edgeIndices)
    {
        for (const Index boundaryEdgeIndex : edgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[boundaryEdgeIndex];
            const std::array<T, N>& edgeNormal = edgeProperty.normal;
            Vec<T, N+2> U(0);
            if (edgeProperty.connectingCells.first >= mesh.numCells)
            {
                // U = mesh.cells.cellAverages[edgeProperty.connectingCells.second];
                U = edgeProperty.reconstructedPushedU.second;
            }
            else if (edgeProperty.connectingCells.second >= mesh.numCells)
            {
                // U = mesh.cells.cellAverages[edgeProperty.connectingCells.first];
                U = edgeProperty.reconstructedPushedU.first;
            }
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
                edgeProperty.numericalFlux = RoeFlux<T, N>::computeNumericalFlux(U, U, directionVectors);
            }
        }
    }
    static void SolidWallConditionAtEdges(Mesh<T, N>& mesh, const std::vector<Index>& edgeIndices)
    {
        for (const Index boundaryEdgeIndex : edgeIndices)
        {
            EdgeProperty<T, N>& edgeProperty = mesh.edges.edgeProperty[boundaryEdgeIndex];
            edgeProperty.numericalFlux = Vec<T, N+2>(0);
            // continue;
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
            if (edgeProperty.connectingCells.first >= mesh.numCells)
            {
                // const Vec<T, N+2>& UR = mesh.cells.cellAverages[edgeProperty.connectingCells.second];
                const Vec<T, N+2>& UR = edgeProperty.reconstructedPushedU.second;
                Vec<T, N> rhoVR(0);
                for (Index i=0; i<N; i++) rhoVR[i] = UR[1+i];
                const T VNormal = Vec<T, N>::dot(rhoVR, edgeNormal);
                // T normalSq = 0;
                // for (Index i=0; i<N; i++) normalSq += edgeNormal[i] * edgeNormal[i];
                // assert(std::abs(normalSq - 1) < EPSILON2);
                Vec<T, N+2> UL = UR;
                for (uint i=0; i<N; i++)
                {
                    UL[1+i] -= 2 * VNormal * edgeNormal[i] ;
                }
                edgeProperty.numericalFlux = RoeFlux<T, N>::computeNumericalFlux(UL, UR, directionVectors);
            }
            else if (edgeProperty.connectingCells.second >= mesh.numCells)
            {
                // const Vec<T, N+2>& UL = mesh.cells.cellAverages[edgeProperty.connectingCells.first];
                const Vec<T, N+2>& UL = edgeProperty.reconstructedPushedU.first;
                Vec<T, N> rhoVL(0);
                for (Index i=0; i<N; i++) rhoVL[i] = UL[1+i];
                const T VNormal = Vec<T, N>::dot(rhoVL, edgeNormal);
                // T normalSq = 0;
                // for (Index i=0; i<N; i++) normalSq += edgeNormal[i] * edgeNormal[i];
                Vec<T, N+2> UR = UL;
                for (uint i=0; i<N; i++)
                {
                    UR[1+i] -= 2 * VNormal * edgeNormal[i] ;
                }
                edgeProperty.numericalFlux = RoeFlux<T, N>::computeNumericalFlux(UL, UR, directionVectors);
                // std::cout<<"Momentum flux: "<<edgeProperty.numericalFlux[N+1]<<std::endl;
                // std::cout<<"Norm of numerical flux: "<<edgeProperty.numericalFlux.norm()<<std::endl;
                // for (uint fieldIdx=0; fieldIdx<N+2; fieldIdx++)
                // {
                //     std::cout<< fieldIdx << "th component of numerical flux: "<<edgeProperty.numericalFlux[fieldIdx]<<std::endl;
                // }
            }

        }
    }
};