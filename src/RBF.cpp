
#include "RBF.hpp"


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
            const Vec<T, N> ubVec = mesh.edges.vertices[cellEdgeIndices[1]][0] - rbfCenter;
            const Vec<T, N> lbVec = mesh.edges.vertices[cellEdgeIndices[0]][0] - rbfCenter;
            T ub = ubVec[0];
            T lb = lbVec[0];
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
            ConnectivityMatrix(rbfId, cellId) = res;
        }
    }
}

template <typename T, uint N>
void RBFs<T, N>::computeRBFCoefficients(const std::vector<T>& cellAveragesX)
{
    assert(N==1);
    assert(cellAveragesX.size() == numCells);
    MatrixSolver<T>::gaussSeidelSolver(ConnectivityMatrix, cellAveragesX, coefficients, 1000, 1e-6);
};

template <typename T, uint N>
T RBFs<T, N>::value(const Vec<T, N>& point)
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