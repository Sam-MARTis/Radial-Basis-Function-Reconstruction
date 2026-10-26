#pragma once

#include "vec.hpp"
#include <algorithm>
#include <array>
#include <cassert>
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

namespace rbfDetail
{
    // 16-point Gauss-Legendre quadrature; exact for the polynomial profiles.
    inline constexpr double gaussNodes[8] = {
        0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
        0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499};
    inline constexpr double gaussWeights[8] = {
        0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
        0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};

    template <typename T, uint N>
    T integrateNumerically(const Function<T, N>& function, T lb, T ub)
    {
        if (lb > ub) std::swap(lb, ub);
        const T radius = function.influenceRegion();
        lb = std::clamp(lb, static_cast<T>(0), radius);
        ub = std::clamp(ub, static_cast<T>(0), radius);
        if (!(ub > lb)) return static_cast<T>(0);

        constexpr int panels = 4;
        const T panelWidth = (ub - lb) / static_cast<T>(panels);
        T sum = static_cast<T>(0);
        for (int panel = 0; panel < panels; ++panel)
        {
            const T lower = lb + panelWidth * static_cast<T>(panel);
            const T upper = lower + panelWidth;
            const T halfWidth = static_cast<T>(0.5) * panelWidth;
            const T midpoint = static_cast<T>(0.5) * (lower + upper);
            for (int i = 0; i < 8; ++i)
            {
                Vec<T, N> plus;
                Vec<T, N> minus;
                plus[0] = midpoint + halfWidth * static_cast<T>(gaussNodes[i]);
                minus[0] = midpoint - halfWidth * static_cast<T>(gaussNodes[i]);
                sum += static_cast<T>(gaussWeights[i]) *
                       (function.value(plus) + function.value(minus));
            }
        }
        return static_cast<T>(0.5) * panelWidth * sum;
    }
}

template<typename T, uint N>
class PolynomialFunction: public Function<T, N>
{
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
    [[deprecated]] T integrate(T, T) const override
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
    [[deprecated]] T integrate(T, T) const override
    {
        return T{};
    }
};

template<typename T, uint N>
class WendlandFunction: public Function<T, N>
{
    const T invInfluenceRadius;
    const T influenceRadius;
public:
    explicit WendlandFunction(const T _influenceRadius)
        : invInfluenceRadius(1.0/_influenceRadius), influenceRadius(_influenceRadius){}
    T value(const Vec<T, N>& p) const override
    {
        const T r = p.norm() * invInfluenceRadius;
        if (r > 1.0) return 0;
        const T a = 1-r;
        return a*a*a*a * (-4*a + 5);
    }
    T influenceRegion() const override
    {
        return influenceRadius;
    }
    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T drNorm = std::max(p.norm(), static_cast<T>(EPSILON1));
        const T r = drNorm * invInfluenceRadius;
        if (r > 1.0) return Vec<T, N>(0);
        const T a = 1-r;
        const T factor = invInfluenceRadius * (-20*a*a*a + 20*a*a*a*a);
        return factor * p/drNorm;
    }
    T integrate(T lb, T ub) const override
    {
        return rbfDetail::integrateNumerically(*this, lb, ub);
    }
};

template <typename T, uint N>
class CompactPolynomialFunction : public Function<T, N>
{
    const T invInfluenceRadius;
    const T influenceRadius;
    const std::vector<T> coefficients;

public:
    CompactPolynomialFunction(const T _influenceRadius, const std::initializer_list<T> _coefficients)
        : invInfluenceRadius(static_cast<T>(1) / _influenceRadius),
          influenceRadius(_influenceRadius),
          coefficients(_coefficients)
    {
        assert(_influenceRadius > static_cast<T>(0));
    }

    T value(const Vec<T, N>& p) const override
    {
        const T r = p.norm() * invInfluenceRadius;
        if (r >= static_cast<T>(1)) return static_cast<T>(0);
        T result = static_cast<T>(0);
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it)
        {
            result = result * (static_cast<T>(1) - r) + *it;
        }
        return result;
    }

    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T distance = p.norm();
        const T r = distance * invInfluenceRadius;
        if (r >= static_cast<T>(1) || distance <= static_cast<T>(EPSILON1))
        {
            return Vec<T, N>(0);
        }
        T derivative = static_cast<T>(0);
        for (std::size_t k = coefficients.size(); k-- > 1;)
        {
            derivative = derivative * (static_cast<T>(1) - r) + static_cast<T>(k) * coefficients[k];
        }
        return p * (-derivative * invInfluenceRadius / distance);
    }

    T influenceRegion() const override
    {
        return influenceRadius;
    }

    T integrate(T lb, T ub) const override
    {
        return rbfDetail::integrateNumerically(*this, lb, ub);
    }
};

template <typename T, uint N>
class WendlandC0Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandC0Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(influenceRadius, {0, 0, 1}) {}
};

template <typename T, uint N>
class WendlandC4Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandC4Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(influenceRadius, {0, 0, 0, 0, 0, 0, 56, -88, 35}) {}
};

template <typename T, uint N>
class WendlandC6Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandC6Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(
              influenceRadius, {0, 0, 0, 0, 0, 0, 0, 0, 66, -154, 121, -32}) {}
};

template <typename T, uint N>
class WendlandD4C4Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WendlandD4C4Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(influenceRadius, {0, 0, 0, 0, 0, 0, 10, -12, 3}) {}
};

template <typename T, uint N>
class WuC2Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WuC2Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(influenceRadius, {0, 0, 0, 0, 32, -40, 12}) {}
};

template <typename T, uint N>
class WuC5Function : public CompactPolynomialFunction<T, N>
{
public:
    explicit WuC5Function(const T influenceRadius)
        : CompactPolynomialFunction<T, N>(
              influenceRadius, {0, 0, 0, 0, 0, 126, -231, 153, -45, 5}) {}
};

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
        if (r >= static_cast<T>(1)) return static_cast<T>(0);
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
        const T derivative = -static_cast<T>(2) * r / (y * y) * profile;
        return p * (derivative * invInfluenceRadius / distance);
    }

    T influenceRegion() const override
    {
        return influenceRadius;
    }

    T integrate(T lb, T ub) const override
    {
        return rbfDetail::integrateNumerically(*this, lb, ub);
    }
};

template <typename T, uint N>
class CompactExponentialFunction : public Function<T, N>
{
    const T invInfluenceRadius;
    const T influenceRadius;
    static constexpr T kappa = static_cast<T>(6);

public:
    explicit CompactExponentialFunction(const T _influenceRadius)
        : invInfluenceRadius(static_cast<T>(1) / _influenceRadius), influenceRadius(_influenceRadius)
    {
        assert(_influenceRadius > static_cast<T>(0));
    }

    T value(const Vec<T, N>& p) const override
    {
        const T r = p.norm() * invInfluenceRadius;
        if (r >= static_cast<T>(1)) return static_cast<T>(0);
        const T norm = static_cast<T>(1) - std::exp(-kappa);
        return (std::exp(-kappa * r) - std::exp(-kappa)) / norm;
    }

    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T distance = p.norm();
        const T r = distance * invInfluenceRadius;
        if (r >= static_cast<T>(1) || distance <= static_cast<T>(EPSILON1))
        {
            return Vec<T, N>(0);
        }
        const T derivative = -kappa * std::exp(-kappa * r) / (static_cast<T>(1) - std::exp(-kappa));
        return p * (derivative * invInfluenceRadius / distance);
    }

    T influenceRegion() const override
    {
        return influenceRadius;
    }

    T integrate(T lb, T ub) const override
    {
        return rbfDetail::integrateNumerically(*this, lb, ub);
    }
};

template <typename T, uint N>
class CompactRationalFunction : public Function<T, N>
{
    const T invInfluenceRadius;
    const T influenceRadius;
    static constexpr T kappa = static_cast<T>(4);

public:
    explicit CompactRationalFunction(const T _influenceRadius)
        : invInfluenceRadius(static_cast<T>(1) / _influenceRadius), influenceRadius(_influenceRadius)
    {
        assert(_influenceRadius > static_cast<T>(0));
    }

    T value(const Vec<T, N>& p) const override
    {
        const T r = p.norm() * invInfluenceRadius;
        if (r >= static_cast<T>(1)) return static_cast<T>(0);
        return (static_cast<T>(1) - r) / (static_cast<T>(1) + kappa * r);
    }

    Vec<T, N> gradient(const Vec<T, N>& p) const override
    {
        const T distance = p.norm();
        const T r = distance * invInfluenceRadius;
        if (r >= static_cast<T>(1) || distance <= static_cast<T>(EPSILON1))
        {
            return Vec<T, N>(0);
        }
        const T denominator = static_cast<T>(1) + kappa * r;
        const T derivative = -(static_cast<T>(1) + kappa) / (denominator * denominator);
        return p * (derivative * invInfluenceRadius / distance);
    }

    T influenceRegion() const override
    {
        return influenceRadius;
    }

    T integrate(T lb, T ub) const override
    {
        return rbfDetail::integrateNumerically(*this, lb, ub);
    }
};

inline const char* rbfKernelName(const RBF_KERNEL_KIND kind)
{
    switch (kind)
    {
        case RBF_KERNEL_KIND::WENDLAND_C0: return "Wendland C0";
        case RBF_KERNEL_KIND::WENDLAND_C2: return "Wendland C2";
        case RBF_KERNEL_KIND::WENDLAND_C4: return "Wendland C4";
        case RBF_KERNEL_KIND::WENDLAND_C6: return "Wendland C6";
        case RBF_KERNEL_KIND::WU_C2: return "Wu C2";
        case RBF_KERNEL_KIND::BUMP_CINF: return "Bump C-inf";
        case RBF_KERNEL_KIND::WENDLAND_D4_C4: return "Wendland D4 C4";
        case RBF_KERNEL_KIND::WU_C5: return "Wu C5";
        case RBF_KERNEL_KIND::COMPACT_EXP: return "Compact exp";
        case RBF_KERNEL_KIND::COMPACT_RATIONAL: return "Compact rational";
        default: return "unknown";
    }
}

inline const char* rbfKernelShortName(const RBF_KERNEL_KIND kind)
{
    switch (kind)
    {
        case RBF_KERNEL_KIND::WENDLAND_C0: return "wendland_c0";
        case RBF_KERNEL_KIND::WENDLAND_C2: return "wendland_c2";
        case RBF_KERNEL_KIND::WENDLAND_C4: return "wendland_c4";
        case RBF_KERNEL_KIND::WENDLAND_C6: return "wendland_c6";
        case RBF_KERNEL_KIND::WU_C2: return "wu_c2";
        case RBF_KERNEL_KIND::BUMP_CINF: return "bump_cinf";
        case RBF_KERNEL_KIND::WENDLAND_D4_C4: return "wendland_d4_c4";
        case RBF_KERNEL_KIND::WU_C5: return "wu_c5";
        case RBF_KERNEL_KIND::COMPACT_EXP: return "compact_exp";
        case RBF_KERNEL_KIND::COMPACT_RATIONAL: return "compact_rational";
        default: return "unknown";
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
        case RBF_KERNEL_KIND::WENDLAND_D4_C4:
            return std::make_unique<WendlandD4C4Function<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::WU_C5:
            return std::make_unique<WuC5Function<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::COMPACT_EXP:
            return std::make_unique<CompactExponentialFunction<T, N>>(influenceRadius);
        case RBF_KERNEL_KIND::COMPACT_RATIONAL:
            return std::make_unique<CompactRationalFunction<T, N>>(influenceRadius);
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
            RBF_KERNEL_KIND::BUMP_CINF,
            RBF_KERNEL_KIND::WENDLAND_D4_C4,
            RBF_KERNEL_KIND::WU_C5,
            RBF_KERNEL_KIND::COMPACT_EXP,
            RBF_KERNEL_KIND::COMPACT_RATIONAL};
}

namespace functionMath
{
    template <typename T, uint N>
    [[deprecated]] T integrateOverTriangle(const Function<T, N>&, const std::array<Vec<T, N>, 3>&, const uint)
    {
        std::cerr << "Warning: integrateOverTriangle is deprecated and not implemented." << std::endl;
        return T{};
    }
}
