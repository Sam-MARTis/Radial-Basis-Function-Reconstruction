#pragma once

#include <omp.h>

template <typename T, uint N>
void RBFs<T, N>::computeConnectivityMatrix(const Mesh<T, N> &mesh)
{
    assert(N==1); // Only implemented for 1D for now
    #pragma omp parallel for
    for (uint rbfId=0; rbfId<numRBFs; rbfId++)
    {

        const Vec<T, N>& rbfCenter = kernelCenters[rbfId];
        for (uint cellId=0; cellId<mesh.numCells; cellId++)
        {

            T res = 0;
            const std::array<uint, static_cast<uint>(1) << N>& cellEdgeIndices = mesh.cells.edgeIndices[cellId];
            T ub = mesh.edges.vertices[cellEdgeIndices[1]][0][0] - rbfCenter[0];
            T lb = mesh.edges.vertices[cellEdgeIndices[0]][0][0] - rbfCenter[0];

            T ubAbs = absVal(ub);
            T lbAbs = absVal(lb);
            if (lb*ub <0) // They have different signs
            {
                ubAbs = ubAbs > kernelFunction->influenceRegion() ? kernelFunction->influenceRegion() : ubAbs;
                lbAbs = lbAbs > kernelFunction->influenceRegion() ? kernelFunction->influenceRegion() : lbAbs;
                res = kernelFunction->integrate(0, ubAbs) + kernelFunction->integrate(0, lbAbs);

            }
            else
            {
                if (lbAbs > ubAbs)
                {
                    std::swap(ub, lb);
                    std::swap(ubAbs, lbAbs);
                }
                if (lbAbs >= kernelFunction->influenceRegion())
                {
                    res = 0;
                }
                else if (ubAbs >= kernelFunction->influenceRegion())
                {
                    res = kernelFunction->integrate(lbAbs, kernelFunction->influenceRegion());
                }
                else
                res = kernelFunction->integrate(lbAbs, ubAbs);
            }
            ConnectivityMatrix(cellId, rbfId) = res/mesh.cells.cellVolumes[cellId];
        }
    }
    if (!MatrixSolver<T>::checkDiagonalDominance(ConnectivityMatrix))
    {
        std::cerr << "Warning: Connectivity matrix is not diagonally dominant. The RBF method may not converge." << std::endl;
    }
}

template <typename T, uint N>
void RBFs<T, N>::computeRBFCoefficients(const std::vector<T>& cellAveragesX, const uint maxIterations, const T tolerance)
{
    assert(N==1);
    // std::cout << "CellAveragesX.size() = " << cellAveragesX.size() << std::endl;
    // std::cout << "numCells = " << numCells << std::endl;
    assert(cellAveragesX.size() == numCells);
    MatrixSolver<T>::gaussSeidelSolver(ConnectivityMatrix, cellAveragesX, coefficients, maxIterations, tolerance);
};

template <typename T, uint N>
void RBFs<T, N>::computeCoefficientDerivatives(const uint maxIterations, const T tolerance)
{
    #pragma omp parallel for
    for (uint cellID=0; cellID<numRBFs; cellID++)
    {
        std::vector<T> b(numCells);
        std::vector<T> x(numRBFs);
        b[cellID] = static_cast<T>(1);
        if (cellID>0) b[cellID-1] = static_cast<T>(0);
        MatrixSolver<T>::gaussSeidelSolver(ConnectivityMatrix, b, x, maxIterations, tolerance);
        std::copy(x.begin(), x.end(), coefficientsDerivatives.data.begin() + cellID*numRBFs);
    }
}

template <typename T, uint N>
void RBFs<T, N>::computeRBFCoefficientsViaDerivatives(const std::vector<T>& cellAveragesX)
{
    for (uint rbfID=0; rbfID<numRBFs; rbfID++)
    {
        T coeff = 0;
        for (uint cellID=0; cellID<numCells; cellID++)
        {
            coeff += coefficientsDerivatives(cellID, rbfID) * cellAveragesX[cellID];
        }
        coefficients[rbfID] = coeff;
    }
}




template <typename T, uint N>
T RBFs<T, N>::value(const Vec<T, N>& point) const
{
    assert(N==1);
    T result = 0;
    for (uint rbfId=0; rbfId<numRBFs; rbfId++)
    {
        const Vec<T, N>& rbfCenter = kernelCenters[rbfId];
        Vec<T, N> diff = point - rbfCenter;
        T dist = absVal(diff[0]);
        if (dist <= kernelFunction->influenceRegion())
        {
            result += coefficients[rbfId] * kernelFunction->value(diff);
        }
    }
    return result;
}
template <typename T, uint N>
T RBFs<T, N>::integrate(T lb, T ub) const
{
    assert(N==1);
    T result = 0;
    for (uint rbfId=0; rbfId<numRBFs; rbfId++)
    {
        T res = 0;
        const T rbfCenter = kernelCenters[rbfId][0];
        const T influenceRadius = kernelFunction->influenceRegion();
        T lbLocal = lb - rbfCenter;
        T ubLocal = ub - rbfCenter;

        T ubAbs = absVal(ubLocal);
        T lbAbs = absVal(lbLocal);
        if (lbLocal*ubLocal <0) // They have different signs
        {
            ubAbs = ubAbs > influenceRadius ? influenceRadius : ubAbs;
            lbAbs = lbAbs > influenceRadius ? influenceRadius : lbAbs;
            res = kernelFunction->integrate(0, ubAbs) + kernelFunction->integrate(0, lbAbs);

        }
        else
        {
            if (lbAbs > ubAbs)
            {
                std::swap(ub, lb);
                std::swap(ubAbs, lbAbs);
            }
            if (lbAbs >= influenceRadius)
            {
                continue;
                res = 0;
            }
            else if (ubAbs >= influenceRadius)
            {
                res = kernelFunction->integrate(lbAbs, influenceRadius);
            }
            else
                res = kernelFunction->integrate(lbAbs, ubAbs);
        }

        result += coefficients[rbfId] * res;
    }
    return result;
}

template <typename T, uint N>
void RBFs<T, N>::computeLocalCellConnectivityMatrix(const Mesh<T, N>& mesh)
{
    assert(N==1); // Only implemented for 1D for now. For 2D -> Gauss integration of triangles
    #pragma omp parallel for
    for (uint topCellID=0; topCellID<numCells; topCellID++)
    {
        Neighbourhood<T, N>& neighbourhood = mesh.cells.neighbourhoods[topCellID];
        const std::vector<Index>& neighbourCellIndices = neighbourhood.neighbourCellIndices;
        const uint numNeighbours = neighbourhood.numNeighbours;

        // The local structural matrix is square: row i is the equation obtained by
        // integrating over neighbour cell i, column j is the RBF centred at neighbour cell j.
        // The neighbourhood includes cellID itself, so the cell that owns this matrix is one
        // of the neighbours and its own basis function is part of the local interpolation.
        Matrix<T> ALocal(numNeighbours, numNeighbours);

        for (uint j=0; j<numNeighbours; j++)
        {
            const Index rbfId = neighbourCellIndices[j];
            const Vec<T, N>& rbfCenter = kernelCenters[rbfId];
            for (uint i=0; i<numNeighbours; i++)
            {
                const Index cellId = neighbourCellIndices[i];
                T res = 0;
                const std::array<uint, static_cast<uint>(1) << N>& cellEdgeIndices = mesh.cells.edgeIndices[cellId];
                T ub = mesh.edges.vertices[cellEdgeIndices[1]][0][0] - rbfCenter[0];
                T lb = mesh.edges.vertices[cellEdgeIndices[0]][0][0] - rbfCenter[0];

                T ubAbs = absVal(ub);
                T lbAbs = absVal(lb);
                if (lb*ub <0) // They have different signs
                {
                    ubAbs = ubAbs > kernelFunction->influenceRegion() ? kernelFunction->influenceRegion() : ubAbs;
                    lbAbs = lbAbs > kernelFunction->influenceRegion() ? kernelFunction->influenceRegion() : lbAbs;
                    res = kernelFunction->integrate(0, ubAbs) + kernelFunction->integrate(0, lbAbs);

                }
                else
                {
                    if (lbAbs > ubAbs)
                    {
                        std::swap(ub, lb);
                        std::swap(ubAbs, lbAbs);
                    }
                    if (lbAbs >= kernelFunction->influenceRegion())
                    {
                        res = 0;
                    }
                    else if (ubAbs >= kernelFunction->influenceRegion())
                    {
                        res = kernelFunction->integrate(lbAbs, kernelFunction->influenceRegion());
                    }
                    else
                        res = kernelFunction->integrate(lbAbs, ubAbs);
                }
                ALocal(i, j) = res/mesh.cells.cellVolumes[cellId];
            }
        }

        neighbourhood.ALocalInv = MatrixSolver<T>::inverse(ALocal);
    }
}

template <typename T, uint N>
void RBFs<T, N>::updateMeshRBFCoefficients(const Mesh<T, N>& mesh)
{
    for (Index cellId = 0; cellId<mesh.numCells; cellId++)
    {
        const Neighbourhood<T, N>& neighbourhood = mesh.cells.neighbourhoods[cellId];
        const std::vector<Index>& neighbourCellIndices = neighbourhood.neighbourCellIndices;
        const uint numNeighbours = neighbourhood.numNeighbours;
        std::vector<Vec<T, N+1>> cellAverages(numNeighbours);
        for (uint neighbourCellIdx=0; neighbourCellIdx<numNeighbours; neighbourCellIdx++)
        {
            const Index neighbourCellId = neighbourCellIndices[neighbourCellIdx];
            cellAverages[neighbourCellIdx] = mesh.cells.cellAverages[neighbourCellId][0];
        }
        // std::array<std::vector<T>, N+2> rbfCoefficients;
        for (uint fieldVarIdx=0; fieldVarIdx<N+2; fieldVarIdx++)
        {
            std::vector<T> coeffs(numNeighbours);
            for (uint neighbourCellIdx=0; neighbourCellIdx<numNeighbours; neighbourCellIdx++)
            {
                coeffs[neighbourCellIdx] = cellAverages[neighbourCellIdx][fieldVarIdx];
            }
            neighbourhood.rbfCoefficients[fieldVarIdx] = MatrixSolver<T>::matrixStdVectorMultiply(neighbourhood.ALocalInv, coeffs);
        }
        
        // mesh.cells.neighbourhoods[cellId].rbfCoefficients = rbfCoefficients;
        // mesh.cells.neighbourhoods[cellId].rbfCoefficients = rbfCoefficients;
    }
}
