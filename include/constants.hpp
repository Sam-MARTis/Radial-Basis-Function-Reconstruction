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
    NO_LIMITER = 0,
    VENKATAKRISHNAN = 1,
    VAN_LEER = 2,
    MINMOD = 3,
    SUPERBEE = 4,
    MC = 5,
    VAN_ALBADA = 6,
    BARTH_JESPERSEN = 7,
    KOREN = 8,
    COUNT = 9,
};

constexpr LIMITER_KIND LIMITER = LIMITER_KIND::VENKATAKRISHNAN;
constexpr bool VAN_LEER_MONOTONE = false;

enum class RBF_KERNEL_KIND : unsigned int
{
    WENDLAND_C0 = 0,
    WENDLAND_C2 = 1,
    WENDLAND_C4 = 2,
    WENDLAND_C6 = 3,
    WU_C2 = 4,
    BUMP_CINF = 5,
    WENDLAND_D4_C4 = 6,
    WU_C5 = 7,
    COMPACT_EXP = 8,
    COMPACT_RATIONAL = 9,
    COUNT = 10,
};

constexpr RBF_KERNEL_KIND DEFAULT_RBF_KERNEL = RBF_KERNEL_KIND::WENDLAND_C4;

constexpr double EIGENVALUE_SMOOTHENING_DELTA = 0.20;
constexpr double SMOOTH_EIGENVALUE_MAGNITUDE(const double lambda)
{
    const double lambdaMagnitude = std::abs(lambda);
    return lambdaMagnitude > EIGENVALUE_SMOOTHENING_DELTA
               ? lambdaMagnitude
               : (1.0 / (2 * EIGENVALUE_SMOOTHENING_DELTA)) *
                     (lambdaMagnitude * lambdaMagnitude + EIGENVALUE_SMOOTHENING_DELTA * EIGENVALUE_SMOOTHENING_DELTA);
}

constexpr double DOMAIN_X_MAX = 10.0;

constexpr uint DEFAULT_NEIGHBOURHOOD_DEPTH = 3;
constexpr uint LIMITER_NEIGHBOURHOOD_DEPTH =
    (LIMITER == LIMITER_KIND::VAN_LEER) ? 1u : DEFAULT_NEIGHBOURHOOD_DEPTH;

constexpr double KERNEL_SUPPORT_MULTIPLIER = 1.3;

constexpr double RBF_FIRST_ORDER_BLEND = 0.00;
constexpr uint USE_RBF_DIRECTLY_FOR_FLUXES = 2;
constexpr double RBF_SECOND_ORDER_BLEND = 1.0;
constexpr double CV = 1.0;
