#pragma once
#include <cassert>

#include "cells.hpp"

template <typename T, uint N>
class Solver
{
    Edges<T, N>& edges;
    Cells<T, N>& cells;

    uint dimension;
    uint numberOfCells;
    uint numberOfEdges;
    T a;
    [[nodiscard]] bool verifyDomainIntegrity() const;
    void calculateFluxes();
    void updateCellAverages(T dt);
    void applyBoundaryConditions();

    public:
    Solver(Edges<T, N>& _edges, Cells<T, N>& _cells, T _a);
    void step(T dt)
    {
        calculateFluxes();
        applyBoundaryConditions();
        updateCellAverages(dt);
    }


    
};

#include "fluid-solver.tpp"
