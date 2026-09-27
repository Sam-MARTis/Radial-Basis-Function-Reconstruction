#pragma once
#include <cstddef>
#include <Eigen/Dense>

using uint = unsigned int;
using Index = uint;
using EigenMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


constexpr double DEFAULT_SCREEN_WIDTH = 800;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;


constexpr double EPSILON1 = 1e-12;
constexpr double EPSILON2 = 1e-8;
constexpr double EPSILON3 = 1e-4;
constexpr double LIMITER_K = 0.001;
constexpr double EIGENVALUE_SMOOTHENING_DELTA = 0.25;
constexpr double SMOOTH_EIGENVALUE_MAGNITUDE(const double lambda)
{
    const double lambdaMagnitude = std::abs(lambda);
    return lambdaMagnitude> EIGENVALUE_SMOOTHENING_DELTA ? lambdaMagnitude : (1.0/(2*EIGENVALUE_SMOOTHENING_DELTA)) * (lambdaMagnitude*lambdaMagnitude + EIGENVALUE_SMOOTHENING_DELTA*EIGENVALUE_SMOOTHENING_DELTA);
}



constexpr double DOMAIN_X_MAX = 10.0;
constexpr uint DEFAULT_NEIGHBOURHOOD_DEPTH = 3;


template <typename T>
constexpr T absVal(T value) {
    return (value < 0) ? -value : value;
}

constexpr double CV = 1.0;