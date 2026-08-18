#include "mesh.hpp"
#include "vec.hpp"
#include "fluid-solver.hpp"

MeshGenerator::MeshGenerator(uint _dimensions, Vec<double, 2> _domainBounds, uint _numberOfCells, MeshType _meshType)
    : dimensions(_dimensions), domainBounds(_domainBounds), numberOfCells(_numberOfCells), meshType(_meshType)
{
    assert(dimensions ==1); // Will implement 2D and 3D later.
    // assert(dimensions > 0 && dimensions <= 3);
    assert(numberOfCells > 0);
    assert(meshType == MeshType::CARTESIAN);
}


Mesh<double, 1> MeshGenerator::generateMesh()
{
    Mesh<double, 1> mesh;
    mesh.numCells = numberOfCells;
    mesh.meshType = meshType;

    double cellSize = (domainBounds[1] - domainBounds[0]) / numberOfCells;
    // double targetPoint = (domainBounds[0] + domainBounds[1]) / 2.0;
    double targetPoint =  3*cellSize;
    for (uint i = 0; i < numberOfCells + 1; i++)
    {
        
        // Create edges
        // if (i < numberOfCells+1)
        // {
            EdgeProperty<double, 1> edge;
            edge.connectingCells = {i==0 ? (0) : i-1, i == numberOfCells ? numberOfCells-1 : i};
            edge.normal = {1.0}; // Normal pointing from left to right
            edge.numericalFlux = 0.0; // Initialize numerical flux to zero
            mesh.edges.edgeProperty.push_back(edge);
            mesh.edges.vertices.push_back({Vec<double, 1>{domainBounds[0] + (i) * cellSize}});
        // }

        if(i==numberOfCells) continue; // Skip the last cell, as it doesn't have a right edge.
        Vec<double, 1> cellCenter{domainBounds[0] + (i + 0.5) * cellSize};
        mesh.cells.cellCentersPositions.push_back(cellCenter);
        mesh.cells.cellVolumes.push_back(cellSize);
        mesh.cells.cellAverages.push_back({static_cast<double>(cellCenter[0]< targetPoint)});
        mesh.cells.edgeIndices.push_back({i, i+1}); 
    }


    mesh.numEdges = mesh.edges.edgeProperty.size();
    return mesh;
}
