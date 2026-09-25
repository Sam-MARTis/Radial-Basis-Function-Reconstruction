#pragma once
#include "constants.hpp"

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <cassert>
#include <array>

template<typename T, uint N>
struct Vec
{
    std::array<T, N> data;
    Vec()
    {
        data.fill(T{});
    }
    explicit Vec(T val)
    {
        data.fill(val);
    }
    Vec(std::initializer_list<T> values)
    {
        // data.fill(T{});
        assert(N == values.size());
        std::copy(values.begin(), values.end(), data.begin());
    }
    T& operator[](uint i)
    {
        return data[i];
    }
    const T& operator[](size_t i) const
    {
        return data[i];
    }
    void operator+=(const Vec& vec)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] += vec.data[i];
        }
    }
    void operator-=(const Vec& vec)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] -= vec.data[i];
        }
    }
    Vec operator+(const Vec& vec) const
    {
        Vec result(*this);
        result += vec;
        return result;
    }
    Vec operator-(const Vec& vec) const
    {
        Vec result(*this);
        result -= vec;
        return result;
    }
    void operator*=(const T val)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] *= val;
        }
    }
    void operator/=(const T val)
    {
        *this *= 1/val;
    }
    Vec operator*(const T val) const
    {
        Vec result(*this);
        result *= val;
        return result;
    }
    friend Vec operator*(T scalar, const Vec& v) {
        return v * scalar;
    }

    Vec operator/(const T val) const
    {
        Vec result(*this);
        result /= val;
        return result;
    }
    [[nodiscard]] T norm() const
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += data[i] * data[i];
        }
        return std::sqrt(sum);
    }
    constexpr static T dot(const Vec& a, const Vec& b)
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += a.data[i] * b.data[i];
        }
        return sum;
    }
    constexpr static T dot(const Vec& a, const std::array<T, N>& b)
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += a.data[i] * b[i];
        }
        return sum;
    }
};

namespace VecMath
{
    template <typename T, uint N, uint M>
    Vec<T, N> matrixMultiply(const Vec<Vec<T, M>, N>& matrix, const Vec<T, M>& vector)
    {
        Vec<T, N> result;
        for (uint i = 0; i < N; ++i)
        {
            result[i] = Vec<T, M>::dot(matrix[i], vector);
        }
        return result;
    }

    template <typename T, uint N>
    Vec<T, N> componentwiseMultiply(const Vec<T, N>& a, const Vec<T, N>& b)
    {
        Vec<T, N> result;
        for (uint i = 0; i < N; ++i)
        {
            result[i] = a[i] * b[i];
        }
        return result;
    }
}



// inline double distanceSquared(const Vec2D& p1, const Vec2D& p2)
// {
//     return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
// }
