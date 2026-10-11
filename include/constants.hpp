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

constexpr LIMITER_KIND LIMITER = LIMITER_KIND::VAN_LEER;

/*
 * Caps van Leer psi at 1. Diffusivity knob only -- it is NOT a boundedness fix.
 */
constexpr bool VAN_LEER_MONOTONE = false;

constexpr double EIGENVALUE_SMOOTHENING_DELTA = 0.20;
constexpr double SMOOTH_EIGENVALUE_MAGNITUDE(const double lambda)
{
    const double lambdaMagnitude = std::abs(lambda);
    return lambdaMagnitude> EIGENVALUE_SMOOTHENING_DELTA ? lambdaMagnitude : (1.0/(2*EIGENVALUE_SMOOTHENING_DELTA)) * (lambdaMagnitude*lambdaMagnitude + EIGENVALUE_SMOOTHENING_DELTA*EIGENVALUE_SMOOTHENING_DELTA);
}



constexpr double DOMAIN_X_MAX = 10.0;
constexpr uint DEFAULT_NEIGHBOURHOOD_DEPTH = 5;
constexpr double KERNEL_SUPPORT_MULTIPLIER = 1.9;


constexpr double RBF_FIRST_ORDER_BLEND = 0.00;

constexpr uint USE_RBF_DIRECTLY_FOR_FLUXES = 2;
constexpr double RBF_SECOND_ORDER_BLEND = 1.5;
constexpr double CV = 1.0;