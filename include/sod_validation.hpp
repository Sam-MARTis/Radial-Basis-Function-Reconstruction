#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "sod_shock_tube.hpp"
#include "sod_simulation.hpp"

namespace sod
{

struct FieldError
{
    std::string name;
    double l1 = 0.0;
    double l2 = 0.0;
    double lInf = 0.0;
    double relL1 = 0.0;
    double bias = 0.0;
    std::vector<double> x;
    std::vector<double> simulated;
    std::vector<double> reference;
    std::vector<double> error;
};

struct Validation
{
    bool valid = false;
    double time = 0.0;
    std::vector<FieldError> fields;
    double massError = 0.0;
    double momentumError = 0.0;
    double energyError = 0.0;
    double massErrorRel = 0.0;
    double momentumErrorRel = 0.0;
    double energyErrorRel = 0.0;
    double minimumDensity = 0.0;
    double minimumPressure = 0.0;
    double maximumMach = 0.0;

    const FieldError* find(const std::string& name) const
    {
        for (const FieldError& field : fields)
        {
            if (field.name == name) return &field;
        }
        return nullptr;
    }
};

class Validator
{
public:
    static Validation validate(const Data& data, const Simulation& simulation)
    {
        Validation result;
        result.time = simulation.time();
        result.minimumDensity = simulation.minimumDensity();
        result.minimumPressure = simulation.minimumPressure();
        result.maximumMach = simulation.maximumMach();

        result.massError = simulation.finalMass() - simulation.initialMass();
        result.momentumError = simulation.finalMomentum() - simulation.initialMomentum();
        result.energyError = simulation.finalEnergy() - simulation.initialEnergy();
        result.massErrorRel = relative(result.massError, simulation.initialMass());
        result.momentumErrorRel = relative(result.momentumError, simulation.initialMomentum());
        result.energyErrorRel = relative(result.energyError, simulation.initialEnergy());

        const struct
        {
            Field field;
            const char* name;
        } fields[] = {
            {Field::DENSITY, "rho"},
            {Field::VELOCITY, "u"},
            {Field::PRESSURE, "p"},
            {Field::MACH, "M"},
        };

        for (const auto& entry : fields)
        {
            FieldError error;
            error.name = entry.name;
            const std::size_t count = data.size();
            error.x = data.x();
            error.simulated.reserve(count);
            error.reference.reserve(count);
            error.error.reserve(count);

            double sumAbs = 0.0;
            double sumSquared = 0.0;
            double sumBias = 0.0;
            double sumReference = 0.0;
            double maxAbs = 0.0;

            for (std::size_t point = 0; point < count; ++point)
            {
                const double reference = data.value(entry.field, point);
                const double simulated = simulation.probe(entry.field, data.x()[point]);
                const double difference = simulated - reference;
                error.simulated.push_back(simulated);
                error.reference.push_back(reference);
                error.error.push_back(difference);
                sumAbs += std::abs(difference);
                sumSquared += difference * difference;
                sumBias += difference;
                sumReference += std::abs(reference);
                maxAbs = std::max(maxAbs, std::abs(difference));
            }

            const double inverseCount = count == 0 ? 0.0 : 1.0 / static_cast<double>(count);
            error.l1 = sumAbs * inverseCount;
            error.l2 = std::sqrt(sumSquared * inverseCount);
            error.lInf = maxAbs;
            error.relL1 = sumReference > 0.0 ? error.l1 / (sumReference * inverseCount) : 0.0;
            error.bias = sumBias * inverseCount;
            result.fields.push_back(std::move(error));
        }

        result.valid = true;
        return result;
    }

private:
    static double relative(const double error, const double reference)
    {
        return std::abs(reference) > EPSILON2 ? error / reference : 0.0;
    }
};

}
