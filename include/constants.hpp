#pragma once
#include <cstddef>
#include <Eigen/Dense>

using uint = unsigned int;
using Index = uint;
using EigenMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


constexpr double DEFAULT_SCREEN_WIDTH = 1600;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;


constexpr double EPSILON1 = 1e-12;
constexpr double EPSILON2 = 1e-8;
constexpr double EPSILON3 = 1e-4;
constexpr double LIMITER_K = 0.0000001;


enum class LIMITER_KIND : unsigned int
{
    NO_LIMITER = 0,       // psi = 1 everywhere (pure high order)
    VENKATAKRISHNAN = 1,  // smooth, epsilon-regularised (original)
    VAN_LEER = 2,         // psi = phi(r)/r with phi(r) = (r + |r|)/(1 + |r|)
};

// Compile-time overrides are used only by the headless comparison harness
// (tools/rbf_kernel_study.cpp) when it is rebuilt with -D... options. The app
// keeps the defaults below.
#ifndef RBF_LIMITER_KIND_OVERRIDE
constexpr LIMITER_KIND LIMITER = LIMITER_KIND::VENKATAKRISHNAN;
#else
constexpr LIMITER_KIND LIMITER = static_cast<LIMITER_KIND>(RBF_LIMITER_KIND_OVERRIDE);
#endif

// Caps van Leer psi at 1. Diffusivity knob, NOT a boundedness fix.
constexpr bool VAN_LEER_MONOTONE = false;

// Compactly supported radial kernels available for the local reconstruction.
// The smoothness label is the continuity of the 1D radial profile at the support
// boundary; higher smoothness gives a more accurate interpolation but a more
// ill-conditioned local system. See include/function.hpp for the definitions.
enum class RBF_KERNEL_KIND : unsigned int
{
    WENDLAND_C0 = 0,   // (1-r)^2
    WENDLAND_C2 = 1,   // (1-r)^4 (1+4r)        -- historical default
    WENDLAND_C4 = 2,   // (1-r)^6 (35r^2+18r+3)
    WENDLAND_C6 = 3,   // (1-r)^8 (32r^3+25r^2+8r+1)
    WU_C2 = 4,         // (1-r)^4 (4+16r+12r^2)
    BUMP_CINF = 5,     // exp(-1/(1-r^2)), C-infinity bump
    COUNT = 6,
};

// Kernel used by src/main.cpp until changed in the UI.
constexpr RBF_KERNEL_KIND DEFAULT_RBF_KERNEL = RBF_KERNEL_KIND::WENDLAND_C4;

constexpr double EIGENVALUE_SMOOTHENING_DELTA = 0.20;
constexpr double SMOOTH_EIGENVALUE_MAGNITUDE(const double lambda)
{
    const double lambdaMagnitude = std::abs(lambda);
    return lambdaMagnitude> EIGENVALUE_SMOOTHENING_DELTA ? lambdaMagnitude : (1.0/(2*EIGENVALUE_SMOOTHENING_DELTA)) * (lambdaMagnitude*lambdaMagnitude + EIGENVALUE_SMOOTHENING_DELTA*EIGENVALUE_SMOOTHENING_DELTA);
}



constexpr double DOMAIN_X_MAX = 10.0;

// RBF fit neighbourhood depth (edge hops); also sets kernel support in src/main.cpp.
#ifndef RBF_NEIGHBOURHOOD_DEPTH_OVERRIDE
constexpr uint DEFAULT_NEIGHBOURHOOD_DEPTH = 10;
#else
constexpr uint DEFAULT_NEIGHBOURHOOD_DEPTH = RBF_NEIGHBOURHOOD_DEPTH_OVERRIDE;
#endif


#ifndef RBF_LIMITER_NEIGHBOURHOOD_DEPTH_OVERRIDE
constexpr uint LIMITER_NEIGHBOURHOOD_DEPTH = (LIMITER == LIMITER_KIND::VAN_LEER) ? 1u : DEFAULT_NEIGHBOURHOOD_DEPTH;
#else
constexpr uint LIMITER_NEIGHBOURHOOD_DEPTH = RBF_LIMITER_NEIGHBOURHOOD_DEPTH_OVERRIDE;
#endif

constexpr double KERNEL_SUPPORT_MULTIPLIER = 1.9;


constexpr double RBF_FIRST_ORDER_BLEND = 0.00; // Zero -> No first order
constexpr uint USE_RBF_DIRECTLY_FOR_FLUXES = 2;
constexpr double RBF_SECOND_ORDER_BLEND = 0.0;
constexpr double CV = 1.0;
