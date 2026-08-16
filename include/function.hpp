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
        return 1.0/invInfluenceRadius;
    }
    T integrate(T lb, T ub) const override
    {
        lb = lb * invInfluenceRadius;
        ub = ub * invInfluenceRadius;
        // ub = ub >1? static_cast<T>(1): ub;
        assert(ub<=1);
        assert(lb <= ub);
        assert(lb >= 0);
        assert(lb <= 1);
        // if (lb > 1) return 0;
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

//
//
// template<typename T, uint N>
// class TargetFunction: public Function<N,T>
// {
//
//     const std::vector<Vec<N,T>> gaussianCenters;
//     const std::vector<T> gaussianSigmas;
//     const std::vector<T> gaussianMags;
//     const std::array<std::vector<T>, N> polynomialCoeffs;
//     const uint numGaussians;
//     std::vector<GaussianFunction<T, N>> gaussians;
//     PolynomialFunction<T, N> polynomial;
//     public:
//     TargetFunction(polynomialCoeffs, gaussianCenters, gaussianSigmas, gaussianMags): polynomialCoeffs(polynomialCoeffs), gaussianCenters(gaussianCenters), gaussianSigmas(gaussianSigmas), gaussianMags(gaussianMags), numGaussians(gaussianCenters.size()){
//         polynomial = PolynomialFunction<T, N>(polynomialCoeffs);
//         for(uint i = 0; i < numGaussians; i++){
//             gaussians.push_back(GaussianFunction<T, N>(gaussianCenters[i], gaussianSigmas[i], gaussianMags[i]));
//         }
//     }
//
//     T value(const Vec<T, N>& p) const override
//     {
//         T result = polynomial.value(p);
//         for(uint i = 0; i < numGaussians; i++){
//             result += gaussians[i].value(p);
//         }
//         return result;
//     }
// };




// template<uint N, typename T = double>
// class InterolatedRBFs: public Function<N,T>
// {
// public:
//     InterolatedRBFs(const std::vector<Vec<N,T>>& points, const std::vector<T>& values, T epsilon = 1.0);
//     T value(const Vec<N,T>& p) const override;
// private:
//     std::vector<Vec<N,T>> points;       