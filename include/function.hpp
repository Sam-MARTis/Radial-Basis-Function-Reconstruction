#pragma once


#include "vec.hpp"
#include <vector>
template<typename T, uint N>
class Function
{
public:
    virtual ~Function() = default;
    virtual T value(const Vec<T, N>& p) const = 0;
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