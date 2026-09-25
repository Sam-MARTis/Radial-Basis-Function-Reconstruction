#pragma once


#include "vec.hpp"
#include <vector>
#include <limits>


template<typename T, uint N>
class Function
{
public:
    virtual ~Function() = default;
    virtual T value(const Vec<T, N>& p) const = 0;
    virtual T influenceRegion() const = 0;
    virtual T integrate(T lb, T ub) const = 0;
};

template<typename T, uint N>
class PolynomialFunction: public Function<T, N>
{
private:
    std::array<std::vector<T>, N> coeffs;
public:
    PolynomialFunction(const std::array<std::vector<T>, N>& _coeffs): coeffs(_coeffs){}
    T value(const Vec<T, N>& p) const override
    {
        T result = 0;
        for(uint i = 0; i < N; i++)
        {
            T term = 1;
            for(uint j = 0; j < coeffs[i].size(); j++)
            {
                term *= std::pow(p[i], j) * coeffs[i][j];
            }
            result *= term;
        }
        return result;  
    }
    T influenceRegion() const override
    {
        return std::numeric_limits<T>::infinity();
    }
    [[deprecated]] T integrate(T lb, T ub) const override
    {
        return T{};
    }
};

template<typename T, uint N>
class GaussianFunction: public Function<T, N>
{
    const Vec<T, N> center;
    const T sigma;
    const T mag;
public:
    GaussianFunction(const Vec<T, N>& center, T sigma, T mag): center(center), sigma(sigma), mag(mag){}
    T value(const Vec<T, N>& p) const override
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += std::pow(p[i] - center[i], 2);
        }
        return mag*std::exp(-sum / (2 * sigma * sigma));
    }
    T influenceRegion() const override
    {
        return std::numeric_limits<T>::infinity();
    }
    [[deprecated]] T integrate(T lb, T ub) const override
    {
        return T{};
    }
};

template<typename T, uint N>
class WendlandFunction: public Function<T, N>
{
    // C^2 wendland function
    const T invInfluenceRadius;
    const T influenceRadius;
public:
    WendlandFunction(const T _influenceRadius): invInfluenceRadius(1.0/_influenceRadius), influenceRadius(_influenceRadius){}
    T value(const Vec<T, N>& p) const override
    {

        const T r = p.norm() * invInfluenceRadius;
        if (r > 1.0) return 0;
        const T a = 1-r;
        return a*a*a*a * (-4*a + 5);
    }
    T influenceRegion() const override
    {
        assert(absVal(invInfluenceRadius - 1.0/influenceRadius) < EPSILON1);
        return influenceRadius;
    }
    T integrate(T lb, T ub) const override
    {
        lb = lb * invInfluenceRadius;
        ub = ub * invInfluenceRadius;
        assert(ub<=1);
        assert(lb <= ub);
        assert(lb >= 0);
        assert(lb <= 1);
        const T y1 = (1-ub);
        const T y0 = (1-lb);
        const T y1pow5 = y1*y1*y1*y1*y1;
        const T y0pow5 = y0*y0*y0*y0*y0;
        const T res1 = 3*y1pow5*(-4*y1 + 5) + 2*y1pow5*y1;
        const T res0 = 3*y0pow5*(-4*y0 + 5) + 2*y0pow5*y0;
        const T res = (-influenceRadius/static_cast<T>(15))*(res1 - res0);
        return res;
    }
};

