#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "mesh.hpp"
#include "fluid-solver.hpp"
#include "RBF.hpp"
#include "sod_shock_tube.hpp"

namespace sod
{

struct Options
{
    uint cells = 200;
    RBF_KERNEL_KIND kernel = RBF_KERNEL_KIND::WENDLAND_C4;
    double supportMultiplier = KERNEL_SUPPORT_MULTIPLIER;
    LIMITER_KIND limiter = LIMITER;
    uint limiterDepth = LIMITER_NEIGHBOURHOOD_DEPTH;
    double cfl = 0.2;
};

class Simulation
{
public:
    Simulation(const Data& data, const Options& options)
        : data_(&data), options_(options)
    {
        const Conditions& conditions = data.conditions();
        if (options_.cells < 4) throw std::runtime_error("at least four cells are required");

        MeshGenerator generator(1, Vec<double, 2>{0.0, conditions.domainLength}, options_.cells,
                                MeshType::CARTESIAN, options_.limiterDepth);
        mesh_ = std::make_unique<Mesh<double, 1>>(generator.generateMesh());
        mesh_->fixEdgesConnectivityOrder();

        HyperCuboid<double, 1> leftBoundary(Vec<double, 1>{0.0}, Vec<double, 1>{0.0});
        HyperCuboid<double, 1> rightBoundary(Vec<double, 1>{conditions.domainLength},
                                             Vec<double, 1>{conditions.domainLength});
        std::vector<std::pair<uint, std::unique_ptr<Domain<double, 1>>>> boundaryDomains;
        boundaryDomains.emplace_back(0, std::make_unique<HyperCuboid<double, 1>>(leftBoundary));
        boundaryDomains.emplace_back(1, std::make_unique<HyperCuboid<double, 1>>(rightBoundary));
        mesh_->identifyBoundaryEdgesAndCells(boundaryDomains);

        const Vec<double, 3> left = data.leftConserved();
        const Vec<double, 3> right = data.rightConserved();
        for (uint cell = 0; cell < options_.cells; ++cell)
        {
            mesh_->cells.cellAverages[cell] = mesh_->cells.cellCentersPositions[cell][0] < conditions.interfaceX
                                                  ? left
                                                  : right;
        }
        initialMass_ = conservedSum(0);
        initialMomentum_ = conservedSum(1);
        initialEnergy_ = conservedSum(2);

        solver_ = std::make_unique<Solver<double, 1>>(*mesh_);
        solver_->setLimiter(options_.limiter);
        solver_->prescribeBoundaryConditions(0, BOUNDARY_CONSTRAINT::EXTRAPOLATED);
        solver_->prescribeBoundaryConditions(1, BOUNDARY_CONSTRAINT::EXTRAPOLATED);

        cellSize_ = conditions.domainLength / static_cast<double>(options_.cells);
        const double supportRadius = options_.supportMultiplier * DEFAULT_NEIGHBOURHOOD_DEPTH * cellSize_;
        rbfs_ = std::make_unique<RBFs<double, 1>>(
            options_.cells, options_.cells, mesh_->cells.cellCentersPositions,
            makeRBFKernel<double, 1>(options_.kernel, supportRadius));
        rbfs_->computeLocalCellConnectivityMatrix(*mesh_);

        const double waveSpeed = std::max(data.maxWaveSpeed(), 1.0);
        dt_ = options_.cfl * cellSize_ / waveSpeed;
        const uint steps = std::max<uint>(1, static_cast<uint>(std::ceil(conditions.time / dt_)));
        dt_ = conditions.time / static_cast<double>(steps);
        steps_ = steps;
    }

    void runToEnd()
    {
        const double endTime = data_->conditions().time;
        for (uint step = 0; step < steps_; ++step)
        {
            solver_->stepRBFBased(*rbfs_, dt_);
            time_ += dt_;
            for (uint cell = 0; cell < mesh_->numCells; ++cell)
            {
                const Vec<double, 3>& state = mesh_->cells.cellAverages[cell];
                if (!std::isfinite(state[0]) || !std::isfinite(state[2]) || state[0] <= 0.0)
                {
                    throw std::runtime_error("Sod simulation diverged");
                }
            }
        }
        time_ = endTime;
    }

    double time() const { return time_; }
    double endTime() const { return data_->conditions().time; }
    double dt() const { return dt_; }
    double cellSize() const { return cellSize_; }
    uint cells() const { return options_.cells; }
    double initialMass() const { return initialMass_; }
    double initialMomentum() const { return initialMomentum_; }
    double initialEnergy() const { return initialEnergy_; }
    double finalMass() const { return conservedSum(0); }
    double finalMomentum() const { return conservedSum(1); }
    double finalEnergy() const { return conservedSum(2); }
    const Mesh<double, 1>& mesh() const { return *mesh_; }

    double pressure(const Vec<double, 3>& state) const
    {
        const double rho = std::max(state[0], EPSILON2);
        const double velocity = state[1] / rho;
        return std::max((data_->conditions().gamma - 1.0) *
                            (state[2] - 0.5 * rho * velocity * velocity),
                        EPSILON2);
    }

    double value(const Field field, const uint cell) const
    {
        const Vec<double, 3>& state = mesh_->cells.cellAverages.at(cell);
        switch (field)
        {
            case Field::DENSITY: return state[0];
            case Field::VELOCITY: return state[1] / std::max(state[0], EPSILON2);
            case Field::PRESSURE: return pressure(state);
            case Field::MACH:
            {
                const double sound = std::sqrt(data_->conditions().gamma * pressure(state) /
                                               std::max(state[0], EPSILON2));
                return std::abs(state[1] / std::max(state[0], EPSILON2)) / std::max(sound, EPSILON2);
            }
        }
        return 0.0;
    }

    double probe(const Field field, const double position) const
    {
        const std::vector<Vec<double, 1>>& centers = mesh_->cells.cellCentersPositions;
        const uint count = mesh_->numCells;
        if (position <= centers.front()[0]) return value(field, 0);
        if (position >= centers.back()[0]) return value(field, count - 1);

        uint low = 0;
        uint high = count - 1;
        while (high - low > 1)
        {
            const uint middle = (low + high) / 2;
            if (centers[middle][0] <= position) low = middle;
            else high = middle;
        }
        const double weight = (position - centers[low][0]) / (centers[high][0] - centers[low][0]);
        return (1.0 - weight) * value(field, low) + weight * value(field, high);
    }

    double minimumDensity() const
    {
        double result = mesh_->cells.cellAverages.front()[0];
        for (const Vec<double, 3>& state : mesh_->cells.cellAverages)
        {
            result = std::min(result, state[0]);
        }
        return result;
    }

    double minimumPressure() const
    {
        double result = pressure(mesh_->cells.cellAverages.front());
        for (const Vec<double, 3>& state : mesh_->cells.cellAverages)
        {
            result = std::min(result, pressure(state));
        }
        return result;
    }

    double maximumMach() const
    {
        double result = 0.0;
        for (uint cell = 0; cell < mesh_->numCells; ++cell)
        {
            result = std::max(result, value(Field::MACH, cell));
        }
        return result;
    }

private:
    double conservedSum(const uint component) const
    {
        double total = 0.0;
        for (uint cell = 0; cell < mesh_->numCells; ++cell)
        {
            total += mesh_->cells.cellAverages[cell][component] * mesh_->cells.cellVolumes[cell];
        }
        return total;
    }

    const Data* data_ = nullptr;
    Options options_;
    std::unique_ptr<Mesh<double, 1>> mesh_;
    std::unique_ptr<Solver<double, 1>> solver_;
    std::unique_ptr<RBFs<double, 1>> rbfs_;
    double cellSize_ = 0.0;
    double dt_ = 0.0;
    double time_ = 0.0;
    uint steps_ = 0;
    double initialMass_ = 0.0;
    double initialMomentum_ = 0.0;
    double initialEnergy_ = 0.0;
};

}
