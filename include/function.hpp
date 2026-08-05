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


// template<uint N, typename T = double>
// class InterolatedRBFs: public Function<N,T>
// {
// public:
//     InterolatedRBFs(const std::vector<Vec<N,T>>& points, const std::vector<T>& values, T epsilon = 1.0);
//     T value(const Vec<N,T>& p) const override;
// private:
//     std::vector<Vec<N,T>> points;       