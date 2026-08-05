#pragma once
#include "constants.hpp"

#include <cmath>
#include <initializer_list>
#include <stdexcept>

template<uint N, typename T = double>
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
        if (values.size() != N)
            throw std::invalid_argument("Incorrect number of elements");
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
    void operator+=(Vec& vec)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] += vec.data[i];
        }
    }
    void operator-=(Vec& vec)
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

};

// inline double distanceSquared(const Vec2D& p1, const Vec2D& p2)
// {
//     return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
// }
