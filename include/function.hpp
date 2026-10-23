#pragma once


#include "vec.hpp"
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>


template<typename T, uint N>
class Function
{
public:
    virtual ~Function() = default;
    virtual T value(const Vec<T, N>& p) const = 0;
    virtual T influenceRegion() const = 0;
    virtual T integrate(T lb, T ub) const = 0;
    virtual Vec<T, N> gradient(const Vec<T, N>& p) const = 0;
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
    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        Vec<T, N> grad(0);
        for(uint i = 0; i < N; i++)
        {
            T term = 0;
            for(uint j = 1; j < coeffs[i].size(); j++)
            {
                term += j * coeffs[i][j] * std::pow(p[i], j - 1);
            }
            grad[i] = term;
        }
        return grad;
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
    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        Vec<T, N> grad(0);
        T val = value(p);
        for (uint i = 0; i < N; i++)
        {
            grad[i] = -val * (p[i] - center[i]) / (sigma * sigma);
        }
        return grad;
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
        assert(std::abs(invInfluenceRadius - 1.0/influenceRadius) < EPSILON1);
        return influenceRadius;
    }
    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T drNorm = std::max(p.norm(), EPSILON1);

        const T r = drNorm * invInfluenceRadius;
        if (r > 1.0) return Vec<T, N>(0);
        const T a = 1-r;
        const T factor = invInfluenceRadius * (-20*a*a*a + 20*a*a*a*a);
        return factor * p/drNorm;
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


/*
 * Generic compactly supported radial polynomial.
 *
 * The profile is written in the variable y = 1 - r (r = |p| / R):
 *
 *     phi(r) = sum_k c[k] * y^k,        0 <= r <= 1,   phi = 0 outside.
 *
 * Every Wendland and Wu kernel in this file has this form, so the value, the
 * radial derivative and the exact antiderivative can all be evaluated from one
 * coefficient list. Keeping the antiderivative analytic matters: the local
 * cell-average operator is built by integrating the kernel over cells, and a
 * quadrature error there is indistinguishable from a modelling error later.
 */
template <typename T, uint N>
class CompactPolynomialFunction : public Function<T, N>
{
    const T invInfluenceRadius;
    const T influenceRadius;
    const std::vector<T> coefficients; // c[k] for y^k

public:
    CompactPolynomialFunction(const T _influenceRadius, const std::initializer_list<T> _coefficients)
        : invInfluenceRadius(static_cast<T>(1) / _influenceRadius),
          influenceRadius(_influenceRadius),
          coefficients(_coefficients)
    {
        assert(_influenceRadius > static_cast<T>(0));
        assert(!coefficients.empty());
    }

    T value(const Vec<T, N>& p) const override
    {
        const T r = p.norm() * invInfluenceRadius;
        if (r >= static_cast<T>(1))
        {
            return static_cast<T>(0);
        }
        return evaluatePolynomial(static_cast<T>(1) - r);
    }

    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T distance = p.norm();
        const T r = distance * invInfluenceRadius;
        if (r >= static_cast<T>(1))
        {
            return Vec<T, N>(0);
        }
        // At the centre of a radial basis function the direction is undefined. Kernels
        // with a non-zero c[1] have a finite radial slope there, but every kernel used
        // for reconstruction has c[1] = 0, so the limit is the zero vector.
        if (distance <= static_cast<T>(EPSILON1))
        {
            return Vec<T, N>(0);
        }
        const T derivative = -evaluateDerivative(static_cast<T>(1) - r); // d(phi)/dr
        return p * (derivative * invInfluenceRadius / distance);
    }

    T influenceRegion() const override
    {
        return influenceRadius;
    }

    /* Exact integral of the physical (un-normalised) profile over [lb, ub]. */
    T integrate(T lb, T ub) const override
    {
        if (lb > ub)
        {
            std::swap(lb, ub);
        }
        lb = std::clamp(lb, static_cast<T>(0), influenceRadius);
        ub = std::clamp(ub, static_cast<T>(0), influenceRadius);
        return influenceRadius * (evaluatePrimitive(ub * invInfluenceRadius) -
                                  evaluatePrimitive(lb * invInfluenceRadius));
    }

    [[nodiscard]] const std::vector<T>& profileCoefficients() const
    {
        return coefficients;
    }

private:
    /* Horner evaluation of sum_k c[k] y^k. */
    T evaluatePolynomial(const T y) const
    {
        T result = static_cast<T>(0);
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it)
        {
            result = result * y + *it;
        }
        return result;
    }

    /* Derivative with respect to y, i.e. d(phi)/dy = -d(phi)/dr. */
    T evaluateDerivative(const T y) const
    {
        T result = static_cast<T>(0);
        for (std::size_t k = coefficients.size(); k-- > 1;)
        {
            result = result * y + static_cast<T>(k) * coefficients[k];
        }
        return result;
    }

    /* Integral_0^s phi(u) du with u and s in normalised radius. */
    T evaluatePrimitive(const T s) const
    {
        T result = static_cast<T>(0);
        for (std::size_t k = 0; k < coefficients.size(); ++k)
        {
            const T exponent = static_cast<T>(k + 1);
            result += coefficients[k] / exponent *
                      (static_cast<T>(1) - std::pow(static_cast<T>(1) - s, exponent));
        }
        return result;
    }
};

/* (1-r)_+^2. The least smooth kernel in the set; useful as a low-order control. */
template <typename T, uint N>
class WendlandC0Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandC0Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(influenceRadius, {static_cast<T>(0), static_cast<T>(0), static_cast<T>(1)})
    {
    }
};

/* (1-r)_+^6 (35 r^2 + 18 r + 3): C^4 Wendland kernel. */
template <typename T, uint N>
class WendlandC4Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandC4Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(
              influenceRadius,
              {static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
               static_cast<T>(0), static_cast<T>(0), static_cast<T>(56), static_cast<T>(-88),
               static_cast<T>(35)})
    {
    }
};

/* (1-r)_+^8 (32 r^3 + 25 r^2 + 8 r + 1): C^6 Wendland kernel. */
template <typename T, uint N>
class WendlandC6Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandC6Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(
              influenceRadius,
              {static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
               static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
               static_cast<T>(66), static_cast<T>(-154), static_cast<T>(121), static_cast<T>(-32)})
    {
    }
};

/* Wu's compactly supported C^2 kernel: (1-r)_+^4 (4 + 16 r + 12 r^2). */
template <typename T, uint N>
class WuC2Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WuC2Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(
              influenceRadius,
              {static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
               static_cast<T>(32), static_cast<T>(-40), static_cast<T>(12)})
    {
    }
};

namespace rbfDetail
{
    template <typename T>
    T bumpProfile(const T s)
    {
        if (s >= static_cast<T>(1))
        {
            return static_cast<T>(0);
        }
        const T y = static_cast<T>(1) - s * s;
        return std::exp(-static_cast<T>(1) / y);
    }

    template <typename T, typename F>
    T simpson(const F& f, const T a, const T b)
    {
        const T midpoint = static_cast<T>(0.5) * (a + b);
        return (b - a) / static_cast<T>(6) *
               (f(a) + static_cast<T>(4) * f(midpoint) + f(b));
    }

    template <typename T, typename F>
    T adaptiveSimpsonRecursive(const F& f, const T a, const T b, const T tolerance,
                               const T whole, const int depth)
    {
        const T midpoint = static_cast<T>(0.5) * (a + b);
        const T left = simpson(f, a, midpoint);
        const T right = simpson(f, midpoint, b);
        const T delta = left + right - whole;
        if (depth <= 0 || std::abs(delta) <= static_cast<T>(15) * tolerance)
        {
            return left + right + delta / static_cast<T>(15);
        }
        return adaptiveSimpsonRecursive(f, a, midpoint, tolerance / static_cast<T>(2), left, depth - 1) +
               adaptiveSimpsonRecursive(f, midpoint, b, tolerance / static_cast<T>(2), right, depth - 1);
    }

    template <typename T, typename F>
    T adaptiveSimpson(const F& f, const T a, const T b, const T tolerance, const int maxDepth)
    {
        if (b <= a)
        {
            return static_cast<T>(0);
        }
        return adaptiveSimpsonRecursive(f, a, b, tolerance, simpson(f, a, b), maxDepth);
    }
}

/*
 * C-infinity compact bump: phi(r) = exp(-1 / (1 - r^2)) for r < 1, 0 outside.
 * It has no elementary antiderivative, so integrate() uses adaptive Simpson.
 * The profile is flat to all orders at the support boundary, which is the
 * opposite end of the smoothness spectrum from the polynomial kernels.
 */
template <typename T, uint N>
class BumpCInfFunction : public Function<T, N>
{
    const T invInfluenceRadius;
    const T influenceRadius;

public:
    explicit BumpCInfFunction(const T _influenceRadius)
        : invInfluenceRadius(static_cast<T>(1) / _influenceRadius), influenceRadius(_influenceRadius)
    {
        assert(_influenceRadius > static_cast<T>(0));
    }

    T value(const Vec<T, N>& p) const override
    {
        const T r = p.norm() * invInfluenceRadius;
        if (r >= static_cast<T>(1))
        {
            return static_cast<T>(0);
        }
        const T y = static_cast<T>(1) - r * r;
        return std::exp(-static_cast<T>(1) / y);
    }

    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T distance = p.norm();
        const T r = distance * invInfluenceRadius;
        if (r >= static_cast<T>(1) || distance <= static_cast<T>(EPSILON1))
        {
            return Vec<T, N>(0);
        }
        const T y = static_cast<T>(1) - r * r;
        const T profile = std::exp(-static_cast<T>(1) / y);
        // d/dr exp(-1/(1-r^2)) = -2 r / (1-r^2)^2 * profile.
        const T derivative = -static_cast<T>(2) * r / (y * y) * profile;
        return p * (derivative * invInfluenceRadius / distance);
    }

    T influenceRegion() const override
    {
        return influenceRadius;
    }

    T integrate(T lb, T ub) const override
    {
        if (lb > ub)
        {
            std::swap(lb, ub);
        }
        lb = std::clamp(lb, static_cast<T>(0), influenceRadius);
        ub = std::clamp(ub, static_cast<T>(0), influenceRadius);
        const T a = lb * invInfluenceRadius;
        const T b = ub * invInfluenceRadius;
        const T normalized = rbfDetail::adaptiveSimpson<T>(
            [](const T s) { return rbfDetail::bumpProfile(s); },
            a, b, static_cast<T>(1e-13), 40);
        return influenceRadius * normalized;
    }
};

/* Name used in reports and in the ImGui kernel selector. */
inline const char* rbfKernelName(const RBF_KERNEL_KIND kind)
{
    switch (kind)
    {
        case RBF_KERNEL_KIND::WENDLAND_C0: return "Wendland C0";
        case RBF_KERNEL_KIND::WENDLAND_C2: return "Wendland C2";
        case RBF_KERNEL_KIND::WENDLAND_C4: return "Wendland C4";
        case RBF_KERNEL_KIND::WENDLAND_C6: return "Wendland C6";
        case RBF_KERNEL_KIND::WU_C2:       return "Wu C2";
        case RBF_KERNEL_KIND::BUMP_CINF:   return "Bump C-inf";
        default:                           return "unknown";
    }
}

/* Lower-case, filename-safe name. */
inline const char* rbfKernelShortName(const RBF_KERNEL_KIND kind)
{
    switch (kind)
    {
        case RBF_KERNEL_KIND::WENDLAND_C0: return "wendland_c0";
        case RBF_KERNEL_KIND::WENDLAND_C2: return "wendland_c2";
        case RBF_KERNEL_KIND::WENDLAND_C4: return "wendland_c4";
        case RBF_KERNEL_KIND::WENDLAND_C6: return "wendland_c6";
        case RBF_KERNEL_KIND::WU_C2:       return "wu_c2";
        case RBF_KERNEL_KIND::BUMP_CINF:   return "bump_cinf";
        default:                           return "unknown";
    }
}

template <typename T, uint N>
std::unique_ptr<Function<T, N>> makeRBFKernel(const RBF_KERNEL_KIND kind, const T influenceRadius)
{
    switch (kind)
    {
        case RBF_KERNEL_KIND::WENDLAND_C0:
            return std::make_unique<WendlandC0Function<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::WENDLAND_C2:
            return std::make_unique<WendlandFunction<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::WENDLAND_C4:
            return std::make_unique<WendlandC4Function<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::WENDLAND_C6:
            return std::make_unique<WendlandC6Function<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::WU_C2:
            return std::make_unique<WuC2Function<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::BUMP_CINF:
            return std::make_unique<BumpCInfFunction<T, N>>(influenceRadius);
        default:
            return std::make_unique<WendlandFunction<T, N>>(influenceRadius);
    }
}

inline std::vector<RBF_KERNEL_KIND> allRBFKernelKinds()
{
    return {RBF_KERNEL_KIND::WENDLAND_C0,
            RBF_KERNEL_KIND::WENDLAND_C2,
            RBF_KERNEL_KIND::WENDLAND_C4,
            RBF_KERNEL_KIND::WENDLAND_C6,
            RBF_KERNEL_KIND::WU_C2,
            RBF_KERNEL_KIND::BUMP_CINF};
}

namespace functionMath
{
    template <typename T, uint N>
    [[deprecated]] T integrateOverTriangle(const Function<T, N>& func, const std::array<Vec<T, N>, 3>& triangleVertices, const uint numGaussOrder)
    {
        std::cerr<< "Warning: integrateOverTriangle is deprecated and not implemented." << std::endl;
        return T{};
    }
}

