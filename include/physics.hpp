#pragma once
#include "vec.hpp"
#include "constants.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

enum class DerivedQuantitiesIndices: uint{
    PRESSURE,
    VELOCITY_MAGNITUDE_SQ,
    ENTHALPY_H,
    COUNT,
};

inline constexpr uint numDerivedQuantities = static_cast<uint>(DerivedQuantitiesIndices::COUNT);



template <typename T, int N>
class EulerEquations
{

public:

    static T getPressureFromEnergyRhoVMagSq(const T E, const T rho, const T VMagSq)
    {
        constexpr T gamma = 1.4;
        return (gamma - 1)*(E - 0.5*rho*VMagSq);
    }
    static T getEFromPressureRhoVMagSq(const T p, const T rho, const T VMagSq)
    {
        constexpr T gamma = 1.4;
        return p/(gamma - 1) + 0.5*rho*VMagSq;
    }
    static T getSoundSpeedFromPressureRho(const T p, const T rho)
    {
        constexpr T gamma = 1.4;
        return std::sqrt(gamma*p/rho);
    }
    static T getPressureFromTemperatureDensity(const T Temp, const T rho)
    {
        constexpr T gamma = 1.4;
        return (gamma-1)*rho*CV*Temp;
    }
    static T getRhoFromPressureTemperature(const T p, const T Temp)
    {
        constexpr T gamma = 1.4;
        return p/((gamma-1)*CV*Temp);
    }


    static Vec<T, N+2> projectToLocalCoordinate(const Vec<T, N+2>& U, const std::array<std::array<T, N>, N>& directionVectors)
    {
        Vec<T, N+2> UProjected(0);
        UProjected[0] = U[0];
        UProjected[N + 1] = U[N + 1];
        for (Index i=0; i<N; i++)
        {
            for (Index j=0; j<N; j++) UProjected[1+i]+=U[1+j]*directionVectors[i][j];
        }
        return UProjected;
    }
    static Vec<T, N+2> projectToGlobalCoordinate(const Vec<T, N+2>& UProjected, const std::array<std::array<T, N>, N>& directionVectors)
    {
        Vec<T, N+2> U(0);
        U[0] = UProjected[0];
        U[N + 1] = UProjected[N + 1];
        for (Index i=0; i<N; i++)
        {
            for (Index j=0; j<N; j++) U[1+i]+=UProjected[1+j]*directionVectors[j][i];
        }
        return U;
    }
    // static Matrix
static Vec<T, N+2> computeFlux(const Vec<T, N+2>& U, std::array<T, numDerivedQuantities>& derivedQuantities)
    {
        constexpr T gamma = 1.4;
        const T rho = U[0];
        const T E = U[N+1];

        assert(rho > 0);
        assert(E > 0);




        // guard against non-positive densities
        // T safeRho = rho <= 0 ? static_cast<T>(1e-8) : rho;
        const T safeRho = std::max(rho, static_cast<T>(EPSILON2));
        const T invRho = static_cast<T>(1.0) / safeRho;
        // const std::array<T, N>& normal = directionVectors[0];
        Vec<T, N> velocity(0);
        for (Index i=0; i<N; i++) velocity[i] = U[1+i] * invRho;


        // for (Index i=0; i<N; i++)
        // {
        //     const std::array<T, N>& direction = directionVectors[i];
        //     projectedVelocities[i] = 0;
        //     for (Index j=0; j<N; j++) projectedVelocities[i] += velocity[j] * direction[j];
        // }


        T VNormal = velocity[0];
        // for (Index i=0; i<N; i++) VNormal+= velocity[i]*normal[i];


        T VMagSq = Vec<T, N>::dot(velocity, velocity);
        // for (Index i=0; i<N; i++) VMagSq+= velocity[i]*velocity[i];

        
        const T p = (gamma - 1)*(E - 0.5*safeRho*VMagSq);
        const T H = (E + p) * invRho;

        Vec<T, N+2> flux(0);
        flux[0] = rho*VNormal;
        for (Index i=0; i<N; i++) flux[1+i] = rho*velocity[i]*VNormal;
        // for (Index i=0; i<N; i++) flux[1+i] = rho*velocity[i]*VNormal + p*normal[i];
        flux[1] += p;
        flux[N+1] =  rho*H*VNormal;
        if (!(rho>0))
        {
            std::cout<<"Warning: non-positive density encountered (rho=" << rho << "). Using safeRho=" << safeRho << " for flux computation. Energy=" << E << ", p=" << p << std::endl;
        }
        // Vec<T, static_cast<uint>(DerivedQuantitiesIndices::COUNT)> derivedQuantities;
        derivedQuantities[static_cast<uint>(DerivedQuantitiesIndices::PRESSURE)] = p;
        derivedQuantities[static_cast<uint>(DerivedQuantitiesIndices::VELOCITY_MAGNITUDE_SQ)] = VMagSq;
        derivedQuantities[static_cast<uint>(DerivedQuantitiesIndices::ENTHALPY_H)] = H;
        return flux;
    }

};


template <typename T, uint N>
class RussanovFlux
{
public:
    static Vec<T, N+2> computeNumericalFlux(const Vec<T, N+2>& UL, const Vec<T, N+2>& UR, const std::array<std::array<T, N>, N>& directionVectors)
    {
        const auto ULProjected = EulerEquations<T, N>::projectToLocalCoordinate(UL, directionVectors);
        const auto URProjected = EulerEquations<T, N>::projectToLocalCoordinate(UR, directionVectors);
        // assert(std::abs(normal.norm() - 1.0) < EPSILON2);
        constexpr T gamma = 1.4;
        std::array<T, numDerivedQuantities> derivedQuantitiesL{};
        std::array<T, numDerivedQuantities> derivedQuantitiesR{};

        
        for (uint i = 0; i < N+2; ++i) {
            if (std::isnan(ULProjected[i]) || std::isinf(ULProjected[i])) {
                std::cerr << "RoeFlux: ULProjected contains NaN/Inf at idx " << i << " value=" << ULProjected[i] << "\n";
            }
            if (std::isnan(URProjected[i]) || std::isinf(URProjected[i])) {
                std::cerr << "RoeFlux: URProjected contains NaN/Inf at idx " << i << " value=" << URProjected[i] << "\n";
            }
        }
        const Vec<T, N+2> FL = EulerEquations<T, N>::computeFlux(ULProjected, derivedQuantitiesL);
        const Vec<T, N+2> FR = EulerEquations<T, N>::computeFlux(URProjected, derivedQuantitiesR);


        T pL = derivedQuantitiesL[static_cast<uint>(DerivedQuantitiesIndices::PRESSURE)];
        T pR = derivedQuantitiesR[static_cast<uint>(DerivedQuantitiesIndices::PRESSURE)];
        T safeRhoL = std::max(ULProjected[0], static_cast<T>(EPSILON2));
        T safeRhoR = std::max(URProjected[0], static_cast<T>(EPSILON2));
        T absVnL = absVal(ULProjected[1] / safeRhoL);

        T absVnR = absVal(URProjected[1] / safeRhoR);

        T cL = std::sqrt(std::max(gamma * pL/safeRhoL, static_cast<T>(EPSILON2)));
        T cR = std::sqrt(std::max(gamma * pR/safeRhoR, static_cast<T>(EPSILON2)));
        const T alpha  = std::max(absVnL + cL, absVnR + cR);
        const Vec<T, N+2> F = 0.5 * (FL + FR) - 0.5 * alpha * (URProjected - ULProjected);
        return EulerEquations<T, N>::projectToGlobalCoordinate(F, directionVectors);
    }
};
//

//
template <typename T, uint N>
class RoeFlux
{
public:
    static Vec<T, N+2> computeNumericalFlux(const Vec<T, N+2>& UL, const Vec<T, N+2>& UR, const std::array<std::array<T, N>, N>& directionVectors)
    {
        constexpr T gamma = 1.4;
        const auto ULProjected = EulerEquations<T, N>::projectToLocalCoordinate(UL, directionVectors);
        const auto URProjected = EulerEquations<T, N>::projectToLocalCoordinate(UR, directionVectors);
        std::array<T, numDerivedQuantities> derivedQuantitiesL{};
        std::array<T, numDerivedQuantities> derivedQuantitiesR{};

        const Vec<T, N+2> FL = EulerEquations<T, N>::computeFlux(ULProjected, derivedQuantitiesL);
        const Vec<T, N+2> FR = EulerEquations<T, N>::computeFlux(URProjected, derivedQuantitiesR);

        const T safeRhoL = std::max(ULProjected[0], static_cast<T>(EPSILON2));
        const T rhoLsqrt = std::sqrt(safeRhoL);
        const T invRhoL = static_cast<T>(1.0)/safeRhoL;
        const T safeRhoR = std::max(URProjected[0], static_cast<T>(EPSILON2));
        const T rhoRsqrt = std::sqrt(safeRhoR);
        const T invRhoR = static_cast<T>(1.0)/safeRhoR;
        const T rhoPAv = rhoLsqrt * rhoRsqrt;
        const T invSumSqrtRho = 1.0/(rhoLsqrt + rhoRsqrt);

        std::array<T, N> projectedPAvVelocities{};
        // Vec<T, N> VL(0);
        // Vec<T, N> VR(0);
        T VMagSqPAv = 0.0;
        const T HPAv_raw = (rhoLsqrt * derivedQuantitiesL[static_cast<uint>(DerivedQuantitiesIndices::ENTHALPY_H)] + rhoRsqrt * derivedQuantitiesR[static_cast<uint>(DerivedQuantitiesIndices::ENTHALPY_H)]) * invSumSqrtRho;
        const T HPAv = std::isnan(HPAv_raw) ? static_cast<T>(0) : HPAv_raw;
        std::array<T, N+2> projectedVelocitiesJump{};
        for (uint i=0; i<N; i++)
        {
            const T VL = ULProjected[1+i] * invRhoL;
            const T VR = URProjected[1+i] * invRhoR;
        // }
        //
        // for (Index i=0; i<N; i++)
        // {
            projectedPAvVelocities[i] = (rhoLsqrt * VL + rhoRsqrt * VR) * invSumSqrtRho;
        // }
        //
        //
        // for (Index i=0; i<N; i++)
        // {
            VMagSqPAv += projectedPAvVelocities[i] * projectedPAvVelocities[i];
        // }
        // for (Index i=0; i<N; i++)
        // {
            projectedVelocitiesJump[i] = VR - VL;
        }
        T acoustArg = (gamma - 1) * (HPAv - 0.5 * VMagSqPAv);
        if (!(acoustArg > 0)) acoustArg = std::max(acoustArg, static_cast<T>(EPSILON2));
        const T aPAv = std::sqrt(acoustArg);
        std::array<Vec<T, N+2>, N+2> eigenVectors{};
        Vec<T, N+2> baseEigenVec(1.0);
        baseEigenVec[N+1] = HPAv - aPAv * projectedPAvVelocities[0];
        for (Index i=0; i<N; i++)
        {
            baseEigenVec[1+i] = projectedPAvVelocities[i];
        }
        eigenVectors[0] = baseEigenVec;
        eigenVectors[0][1] -= aPAv;

        eigenVectors[1] = baseEigenVec;
        eigenVectors[1][N+1] = 0.5 * VMagSqPAv;

        eigenVectors[N+1] = baseEigenVec;
        eigenVectors[N+1][1] += aPAv;
        eigenVectors[N+1][N+1] += 2*aPAv * projectedPAvVelocities[0];

        for (Index i=0; i<N-1; i++)
        {
            eigenVectors[2+i] = Vec<T, N+2>(0);
            eigenVectors[2+i][2+i] = 1.0;
            eigenVectors[2+i][N+1] = projectedPAvVelocities[1+i];
        }
        std::array<T, N+2> eigenValuesAbs{};
        eigenValuesAbs.fill(projectedPAvVelocities[0]); // If there's a bug check this. Lets hope i'm using this properly
        eigenValuesAbs[0] = projectedPAvVelocities[0] - aPAv;
        eigenValuesAbs[N+1] = projectedPAvVelocities[0] + aPAv;
        for (uint i=0; i<N+2; i++)
        {
            eigenValuesAbs[i] = SMOOTH_EIGENVALUE_MAGNITUDE(eigenValuesAbs[i]);
        }
        const T invAPavSq = 1.0 / (aPAv * aPAv);
        const T deltaP = derivedQuantitiesR[static_cast<uint>(DerivedQuantitiesIndices::PRESSURE)] - derivedQuantitiesL[static_cast<uint>(DerivedQuantitiesIndices::PRESSURE)];
        const T deltaRho = URProjected[0] - ULProjected[0];
        std::array<T, N+2> waveStrengths{};
        waveStrengths[0] = 0.5 * invAPavSq * (deltaP - rhoPAv * aPAv * projectedVelocitiesJump[0]);
        waveStrengths[1] = deltaRho - deltaP * invAPavSq;
        for (Index i=0; i<N-1; i++)
        {
            waveStrengths[2+i] = rhoPAv * projectedVelocitiesJump[1+i];
        }
        waveStrengths[N+1] = 0.5 * invAPavSq * (deltaP + rhoPAv * aPAv * projectedVelocitiesJump[0]);

        Vec<T, N+2> flux = FL + FR;
        for (Index i=0; i<N+2; i++)
        {
            flux -= eigenValuesAbs[i] *waveStrengths[i] * eigenVectors[i];
        }
        flux *= 0.5;

        // Project momentum flux back to Cartesian coordinates
        // Vec<T, N+2> fluxCopy = flux;
        //
        // for (Index i = 0; i < N; i++)
        // {
        //     flux[1+i] = 0;
        //
        //     for (Index j = 0; j < N; j++)
        //     {
        //         flux[1+i] +=
        //             fluxCopy[1+j] * directionVectors[j][i];
        //     }
        // }

        return EulerEquations<T, N>::projectToGlobalCoordinate(flux, directionVectors);
    }

};


template <typename T, uint N>
class Limiters
{
public:
    static Vec<T, N> Venkatakrishnan(const Vec<T, N>& Uj, const Vec<T, N>& Ui, const Vec<T, N>& Umin, const Vec<T, N>& Umax, T epsilonSq)
    {
        Vec<T, N> psi(0);
        for (Index fieldIdx = 0; fieldIdx < N; fieldIdx++)
        {
            const T delta = Ui[fieldIdx] - Uj[fieldIdx];
            if (std::abs(delta) < EPSILON2)
            {
                psi[fieldIdx] = 1;
                continue;
            }

            const T deltaP = (delta > 0) ? (Umax[fieldIdx] - Uj[fieldIdx]) : (Umin[fieldIdx] - Uj[fieldIdx]);
            const T deltaM = delta;
            const T num = (deltaP * deltaP + epsilonSq) * deltaM + 2.0 * deltaM * deltaM * deltaP;
            const T den = deltaM * (deltaP * deltaP + 2.0 * deltaM * deltaM + deltaM * deltaP + epsilonSq);
            psi[fieldIdx] = (std::abs(den) > EPSILON2) ? (num / den) : 1.0;
        }
        return psi;
    }

};

enum class ReportedDerivedQuantitiesIndices: uint
{
    PRESSURE,
    VELOCITY_MAGNITUDE_SQ,
    SOUND_SPEED,
    MACH_NUMBER,
    ENTROPY,
    COUNT,
};
template <typename T, uint N>
class DerivedQuantities
{
    Mesh<T, N>& mesh;
    T s_inlet = 0;

public:
    std::vector<std::array<T, static_cast<uint>(ReportedDerivedQuantitiesIndices::COUNT)>> derivedQuantitiesR;
    DerivedQuantities(Mesh<T, N>& mesh, T s_inlet) : mesh(mesh), s_inlet(s_inlet)
    {
        derivedQuantitiesR.resize(mesh.numCells);
    }
    void syncDerivedQuantities()
    {
        constexpr T gamma = 1.4;
        const T safeReferenceEntropy = std::max(std::abs(s_inlet), static_cast<T>(EPSILON2));
        for (Index cellIdx = 0; cellIdx < mesh.numCells; cellIdx++)
        {
            std::array<T, static_cast<uint>(ReportedDerivedQuantitiesIndices::COUNT)> vals{};
            const Vec<T, N+2>& U = mesh.cells.cellAverages[cellIdx];

            const T rho = U[0];
            const T totalEnergy = U[N + 1];
            const T safeRho = std::max(rho, static_cast<T>(EPSILON2));

            T momentumMagnitudeSq = 0;
            for (Index i = 0; i < N; ++i)
            {
                const T momentumComponent = U[1 + i];
                momentumMagnitudeSq += momentumComponent * momentumComponent;
            }

            const T VMagSq = momentumMagnitudeSq / (safeRho * safeRho);
            const T p = EulerEquations<T, N>::getPressureFromEnergyRhoVMagSq(totalEnergy, safeRho, VMagSq);
            const T safePressure = std::max(p, static_cast<T>(EPSILON2));
            const T a = std::sqrt(std::max(gamma * safePressure / safeRho, static_cast<T>(EPSILON2)));
            const T M = a > static_cast<T>(EPSILON2) ? std::sqrt(VMagSq) / a : 0;
            if (M>100) std::cout<<" Mach number for cell "<<cellIdx<<" is "<<M<<std::endl;

            const T rhoTerm = std::pow(safeRho, gamma);
            const T entropyValue = std::log(std::max(safePressure / rhoTerm, static_cast<T>(EPSILON2)));
            const T entropyDelta = entropyValue - std::log(safeReferenceEntropy);

            vals[static_cast<uint>(ReportedDerivedQuantitiesIndices::PRESSURE)] = p;
            vals[static_cast<uint>(ReportedDerivedQuantitiesIndices::VELOCITY_MAGNITUDE_SQ)] = VMagSq;
            vals[static_cast<uint>(ReportedDerivedQuantitiesIndices::SOUND_SPEED)] = a;
            vals[static_cast<uint>(ReportedDerivedQuantitiesIndices::MACH_NUMBER)] = M;
            vals[static_cast<uint>(ReportedDerivedQuantitiesIndices::ENTROPY)] = entropyDelta;
            derivedQuantitiesR[cellIdx] = vals;
        }
    }

};

// #include "mesh.tpp"

