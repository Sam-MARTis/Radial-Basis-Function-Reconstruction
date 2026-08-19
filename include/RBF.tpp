#pragma once


template <typename T, uint N>
void RBFs<T, N>::computeConnectivityMatrix(const Mesh<T, N> &mesh)
{
    assert(N==1); // Only implemented for 1D for now
    for (uint rbfId=0; rbfId<numRBFs; rbfId++)
    {

        const Vec<T, N>& rbfCenter = kernelCenters[rbfId];
        for (uint cellId=0; cellId<mesh.numCells; cellId++)
        {

            T res = 0;
            const std::array<uint, static_cast<uint>(1) << N>& cellEdgeIndices = mesh.cells.edgeIndices[cellId];
            T ub = mesh.edges.vertices[cellEdgeIndices[1]][0][0] - rbfCenter[0];
            T lb = mesh.edges.vertices[cellEdgeIndices[0]][0][0] - rbfCenter[0];
            if (cellId == 30)
            {
                std::cout << ub << " " << lb << std::endl;
            }
            // T ub = ubVec[0];
            // T lb = lbVec[0];
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
            if (cellId == 30) std::cout << res << std::endl;
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
    std::vector<T> b(numCells);
    std::vector<T> x(numRBFs);
    for (uint cellID=0; cellID<numRBFs; cellID++)
    {
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