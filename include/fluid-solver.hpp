#pragma once
#include <cassert>

#include "cells.hpp"
#include "constraints.hpp"
#include "RBF.hpp"
template <typename T, uint N>
class Solver
{
    Mesh<T, N>& mesh;

    uint dimension;
    uint numberOfCells;
    uint numberOfEdges;
    std::vector<std::pair<BOUNDARY_CONSTRAINT, Vec<T, N+2>>> boundaryConditions;
    T a;
    [[nodiscard]] [[deprecated]] bool verifyDomainIntegrity() const;
    void calculateFluxes();
    [[deprecated]] void calculateRBFBasedFluxes(const RBFs<T, N>& rbfs, T dt);
    void updateCellAverages(T dt);
    void applyBoundaryConditions();
    void reconstructSolutionAtEdges(const RBFs<T, N>& rbfs);

    public:
    void updateNeighbourhoodRBFCoeffsAndExtremums();
    void prescribeBoundaryConditions(uint boundaryIndex, BOUNDARY_CONSTRAINT boundaryType, Vec<T, N+2> boundaryValue);
    void prescribeBoundaryConditions(uint boundaryIndex, BOUNDARY_CONSTRAINT boundaryType);
    // void prescribeBoundaryConditions(uint boundaryIndex, BOUNDARY_CONSTRAINT boundaryType);

    Solver(Edges<T, N>& _edges, Cells<T, N>& _cells, T _a);
    void stepRBFBased(const RBFs<T, N>& rbfs, T dt)
    {
        // T dt = 0.1;
        calculateRBFBasedFluxes(rbfs, dt);
        applyBoundaryConditions();
        updateCellAverages(1.0);
    }
    void step(T dt)
    {
        calculateFluxes();
        applyBoundaryConditions();
        updateCellAverages(dt);
    }


    
};

#include "fluid-solver.tpp"
