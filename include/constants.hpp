#pragma once
#include <cstddef>


using uint = unsigned int;

constexpr double DEFAULT_SCREEN_WIDTH = 800;
constexpr double DEFAULT_SCREEN_HEIGHT = 800;




template <typename T>
constexpr T absVal(T value) {
    return (value < 0) ? -value : value;
}