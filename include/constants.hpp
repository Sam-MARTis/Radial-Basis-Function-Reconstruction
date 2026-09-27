#pragma once
#include <cstddef>
#include <Eigen/Dense>

using uint = unsigned int;
using Index = uint;
using EigenMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

constexpr double EPSILON1 = 1e-8;
constexpr double DEFAULT_SCREEN_WIDTH = 800;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;

constexpr double DOMAIN_X_MAX = 10.0;
constexpr uint DEFAULT_NEIGHBOURHOOD_DEPTH = 3;


template <typename T>
constexpr T absVal(T value) {
    return (value < 0) ? -value : value;
}