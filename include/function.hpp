#pragma once
#include "vec.hpp"

template<uint N, typename T>
class Function
{
public:
    virtual ~Function() = default;
    virtual T value(const Vec<N,T>& p) const = 0;
};

template<uint N, typename T>
class PolynomialFunction: public Function<N,T>
{
private:
    std::array<std::vector<T>, N> coeffs;
public:
    PolynomialFunction(const std::array<std::vector<T>, N>& coeffs): coeffs(coeffs){}
    T value(const Vec<N,T>& p) const override
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

template<uint N, typename T>
class GaussianFunction: public Function<N,T>
{
    const Vec<N,T> center;
    const T sigma;
    const T mag;
public:
    GaussianFunction(const Vec<N,T>& center, T sigma, T mag): center(center), sigma(sigma){}
    T value(const Vec<N,T>& p) const override
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += std::pow(p[i] - center[i], 2);
        }
        return mag*std::exp(-sum / (2 * sigma * sigma));
    }
};


template<uint N, typename T>
class TargetFunction: public Function<N,T>
{

    const std::vector<Vec<N,T>> gaussianCenters;
    const std::vector<T> gaussianSigmas;
    const std::vector<T> gaussianMags;
    const std::array<std::vector<T>, N> polynomialCoeffs;
    const uint numGaussians;
    std::vector<GaussianFunction<N,T>> gaussians;
    PolynomialFunction<N,T> polynomial;
    public: 
    TargetFunction(polynomialCoeffs, gaussianCenters, gaussianSigmas, gaussianMags): polynomialCoeffs(polynomialCoeffs), gaussianCenters(gaussianCenters), gaussianSigmas(gaussianSigmas), gaussianMags(gaussianMags), numGaussians(gaussianCenters.size()){
        polynomial = PolynomialFunction<N,T>(polynomialCoeffs);
        for(uint i = 0; i < numGaussians; i++){
            gaussians.push_back(GaussianFunction<N,T>(gaussianCenters[i], gaussianSigmas[i], gaussianMags[i]));
        }   
    }
    
    T value(const Vec<N,T>& p) const override
    {
        T result = polynomial.value(p);
        for(uint i = 0; i < numGaussians; i++){
            result += gaussians[i].value(p);
        }
        return result;
    }
};


// template<uint N, typename T = double>
// class InterolatedRBFs: public Function<N,T>
// {
// public:
//     InterolatedRBFs(const std::vector<Vec<N,T>>& points, const std::vector<T>& values, T epsilon = 1.0);
//     T value(const Vec<N,T>& p) const override;
// private:
//     std::vector<Vec<N,T>> points;       