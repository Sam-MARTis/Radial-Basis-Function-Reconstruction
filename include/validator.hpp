/* Exact-solution validation for the 1D Euler Riemann problem. See developerReference.md. */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "constants.hpp"
#include "mesh.hpp"
#include "vec.hpp"

namespace riemann
{

    
/* Storage layout of the data being validated. */
enum class VariableFormulation
{
    CONSERVATIVE,   // U = {rho, rho*u, E}
    PRIMITIVE,      // W = {rho, u, p}
};

/* Index of each field inside a Vec<T, N+2> state. N == 1 for the supported case. */
enum StateIndex : uint
{
    RHO = 0,        // density
    MOMENTUM = 1,   // rho*u (conserved) or u (primitive)
    ENERGY = 2,     // total energy (conserved) or pressure (primitive)
};

template <typename T>
inline T safeRho(const T rho)
{
    return std::max(rho, static_cast<T>(EPSILON2));
}

/* Conserved vector -> primitive vector {rho, u, p}. */
template <typename T, uint N>
Vec<T, N + 2> toPrimitive(const Vec<T, N + 2>& U)
{
    constexpr T gamma = static_cast<T>(1.4);
    Vec<T, N + 2> W(0);
    const T rho = safeRho(U[RHO]);
    W[RHO] = rho;
    T vMagSq = 0;
    for (uint i = 0; i < N; ++i)
    {
        W[1 + i] = U[1 + i] / rho;
        vMagSq += W[1 + i] * W[1 + i];
    }
    W[N + 1] = (gamma - static_cast<T>(1)) * (U[N + 1] - static_cast<T>(0.5) * rho * vMagSq);
    return W;
}

/* Primitive vector -> conserved vector {rho, rho*u, E}. */
template <typename T, uint N>
Vec<T, N + 2> toConserved(const Vec<T, N + 2>& W)
{
    constexpr T gamma = static_cast<T>(1.4);
    Vec<T, N + 2> U(0);
    const T rho = safeRho(W[RHO]);
    T vMagSq = 0;
    U[RHO] = rho;
    for (uint i = 0; i < N; ++i)
    {
        U[1 + i] = rho * W[1 + i];
        vMagSq += W[1 + i] * W[1 + i];
    }
    U[N + 1] = W[N + 1] / (gamma - static_cast<T>(1)) + static_cast<T>(0.5) * rho * vMagSq;
    return U;
}

/* Normalise a caller-supplied state according to its storage layout. */
template <typename T, uint N>
Vec<T, N + 2> asConserved(const Vec<T, N + 2>& state, const VariableFormulation formulation)
{
    return (formulation == VariableFormulation::CONSERVATIVE) ? state : toConserved<T, N>(state);
}

template <typename T, uint N>
Vec<T, N + 2> asPrimitive(const Vec<T, N + 2>& state, const VariableFormulation formulation)
{
    return (formulation == VariableFormulation::CONSERVATIVE) ? toPrimitive<T, N>(state) : state;
}

// ---------------------------------------------------------------------------
// Field descriptors (drive both the report and the ImPlot series)
// ---------------------------------------------------------------------------

enum class FieldVariable
{
    DENSITY,
    VELOCITY,
    PRESSURE,
    MOMENTUM,
    ENERGY,
    SOUND_SPEED,
    MACH_NUMBER,
    ENTROPY,
};

/* Extract one scalar field from a normalised primitive and conserved state pair. */
template <typename T, uint N>
inline T extractField(const FieldVariable field, const Vec<T, N + 2>& W, const Vec<T, N + 2>& U)
{
    constexpr T gamma = static_cast<T>(1.4);
    const T rho = safeRho(W[RHO]);
    switch (field)
    {
        case FieldVariable::DENSITY:     return W[RHO];
        case FieldVariable::VELOCITY:    return W[1];
        case FieldVariable::PRESSURE:    return W[N + 1];
        case FieldVariable::MOMENTUM:    return U[1];
        case FieldVariable::ENERGY:      return U[N + 1];
        case FieldVariable::SOUND_SPEED: return std::sqrt(std::max(gamma * W[N + 1] / rho, static_cast<T>(EPSILON2)));
        case FieldVariable::MACH_NUMBER:
        {
            const T a = std::sqrt(std::max(gamma * W[N + 1] / rho, static_cast<T>(EPSILON2)));
            return std::abs(W[1]) / a;
        }
        case FieldVariable::ENTROPY:     return W[N + 1] / std::pow(rho, gamma);
    }
    return static_cast<T>(0);
}

struct FieldSpec
{
    FieldVariable variable;
    const char* name;      // short name for ImGui tables / ImPlot labels
    const char* unit;
    bool conserved;        // true if this is one of the conserved variables
};

/* The fields reported for a given formulation. Keep in sync with extractField. */
inline std::vector<FieldSpec> defaultFields(const VariableFormulation formulation)
{
    std::vector<FieldSpec> fields;
    if (formulation == VariableFormulation::PRIMITIVE)
    {
        fields.push_back({FieldVariable::DENSITY,     "rho", "kg/m^3", false});
        fields.push_back({FieldVariable::VELOCITY,    "u",   "m/s",    false});
        fields.push_back({FieldVariable::PRESSURE,    "p",   "Pa",     false});
        fields.push_back({FieldVariable::SOUND_SPEED, "a",   "m/s",    false});
        fields.push_back({FieldVariable::MACH_NUMBER, "M",   "-",      false});
    }
    else
    {
        fields.push_back({FieldVariable::DENSITY,  "rho",    "kg/m^3", true});
        fields.push_back({FieldVariable::MOMENTUM, "rho*u",  "kg/m^2s", true});
        fields.push_back({FieldVariable::ENERGY,   "E",      "J/m^3",  true});
        fields.push_back({FieldVariable::PRESSURE, "p",      "Pa",     false});
        fields.push_back({FieldVariable::VELOCITY, "u",      "m/s",    false});
    }
    return fields;
}





template <typename T>
struct FieldStats
{
    std::string name;
    std::string unit;
    bool conserved = false;

    /* Pointwise comparison at the cell centres. */
    T l1 = 0;            // sum |err| / N
    T l2 = 0;            // sqrt( sum err^2 / N )
    T lInf = 0;          // max |err|
    T relL1 = 0;         // l1 / mean |exact|
    T relL2 = 0;         // l2 / sqrt(mean exact^2)
    T relLInf = 0;       // lInf / max |exact|
    T errorIntegral = 0; // integral( sim - exact ) dx, signed
    T exactIntegral = 0; // integral( exact ) dx, Simpson quadrature

    /* Plot series, all over the cell centres. */
    std::vector<T> x;
    std::vector<T> exact;
    std::vector<T> simulated;
    std::vector<T> error; // simulated - exact
};

template <typename T>
struct ConservationStats
{
    T massError = 0;
    T momentumError = 0;
    T energyError = 0;
    T massErrorRel = 0;
    T momentumErrorRel = 0;
    T energyErrorRel = 0;
    T massSimulated = 0, massExact = 0;
    T momentumSimulated = 0, momentumExact = 0;
    T energySimulated = 0, energyExact = 0;
};

/* Quantities describing the exact solution itself. */
template <typename T>
struct ExactSolutionInfo
{
    bool valid = false;
    T pStar = 0;
    T uStar = 0;
    T rhoStarL = 0;
    T rhoStarR = 0;
    T shockSpeed = 0;      // velocity of the right shock
    T contactSpeed = 0;    // == uStar
    T headSpeed = 0;       // left rarefaction head (leftmost wave)
    T tailSpeed = 0;       // left rarefaction tail
};

/* Everything the UI needs: global scalars, per-field stats and plot series. */
template <typename T>
struct VolumetricData
{
    T time = 0;
    uint numCells = 0;
    VariableFormulation formulation = VariableFormulation::CONSERVATIVE;

    std::vector<T> cellCenters;
    std::vector<FieldStats<T>> fields;
    ConservationStats<T> conservation;
    ExactSolutionInfo<T> exactInfo;

    const FieldStats<T>* find(const std::string& name) const
    {
        for (const FieldStats<T>& f : fields)
        {
            if (f.name == name) return &f;
        }
        return nullptr;
    }
};




template <typename T>
class ExactRiemannSolver
{
public:
    /* Primitive left/right states, the interface position and the gas constant. */
    ExactRiemannSolver(const Vec<T, 3>& leftPrimitive,
                       const Vec<T, 3>& rightPrimitive,
                       const T interfacePosition = static_cast<T>(0.5),
                       const T gamma = static_cast<T>(1.4))
        : gamma_(gamma), x0_(interfacePosition), left_(leftPrimitive), right_(rightPrimitive)
    {
        solve();
    }

    /* Convenience overload taking conserved states. */
    ExactRiemannSolver(const Vec<T, 3>& leftState, const Vec<T, 3>& rightState,
                       const VariableFormulation formulation, const T interfacePosition,
                       const T gamma = static_cast<T>(1.4))
        : ExactRiemannSolver(asPrimitive<T, 1>(leftState, formulation),
                             asPrimitive<T, 1>(rightState, formulation),
                             interfacePosition, gamma)
    {
    }

    [[nodiscard]] bool converged() const { return converged_; }
    [[nodiscard]] T interfacePosition() const { return x0_; }
    [[nodiscard]] const Vec<T, 3>& leftState() const { return left_; }
    [[nodiscard]] const Vec<T, 3>& rightState() const { return right_; }
    [[nodiscard]] T soundSpeedStarLeft() const
    {
        return soundSpeed(left_) * std::pow(std::max(pStar_ / left_[2], static_cast<T>(EPSILON2)),
                                            (gamma_ - static_cast<T>(1)) / (static_cast<T>(2) * gamma_));
    }
    [[nodiscard]] T soundSpeedStarRight() const
    {
        return soundSpeed(right_) * std::pow(std::max(pStar_ / right_[2], static_cast<T>(EPSILON2)),
                                             (gamma_ - static_cast<T>(1)) / (static_cast<T>(2) * gamma_));
    }
    [[nodiscard]] T starPressure() const { return pStar_; }
    [[nodiscard]] T starVelocity() const { return uStar_; }
    [[nodiscard]] T starDensityLeft() const { return rhoStarL_; }
    [[nodiscard]] T starDensityRight() const { return rhoStarR_; }

    /* Exact primitive state at (x, t); initial data for t <= 0. */
    [[nodiscard]] Vec<T, 3> samplePrimitive(const T x, const T t) const
    {
        if (!converged_ || t <= static_cast<T>(0))
        {
            return (x < x0_) ? left_ : right_;
        }
        const T S = (x - x0_) / t;
        const T SContact = uStar_;

        if (S <= SContact)
        {
            // Left wave: a rarefaction when p* < p_L, a shock when p* > p_L.
            if (pStar_ < left_[2])
            {
                // Left rarefaction.
                const T aL = soundSpeed(left_);
                if (S <= left_[1] - aL)
                {
                    return left_;   // untouched left state
                }
                const T SHead = left_[1] - aL;
                const T STail = uStar_ - soundSpeedStarLeft();
                if (S <= SHead || S >= STail)
                {
                    // S <= SHead is the untouched state (already handled above)
                    // S >= STail branch is the uniform star state.
                    return (S <= SHead) ? left_ : starState(left_, rhoStarL_, uStar_);
                }
                // Fan: the left Riemann invariant is constant.
                const T u = static_cast<T>(2) / (gamma_ + static_cast<T>(1)) *
                            (aL + static_cast<T>(0.5) * (gamma_ - static_cast<T>(1)) * left_[1] + S);
                const T a = static_cast<T>(2) / (gamma_ + static_cast<T>(1)) *
                            (aL + static_cast<T>(0.5) * (gamma_ - static_cast<T>(1)) * (left_[1] - S));
                const T rho = left_[0] * std::pow(a / aL, static_cast<T>(2) / (gamma_ - static_cast<T>(1)));
                return Vec<T, 3>{rho, u, left_[2] * std::pow(a / aL, static_cast<T>(2) * gamma_ / (gamma_ - static_cast<T>(1)))};
            }
            return starState(left_, rhoStarL_, uStar_); // left shock
        }

        if (pStar_ > right_[2])
        {
            // Right shock: the star state holds between the contact and the shock and the undisturbed right state is beyond the shock.
            const T shockSpeed = this->shockSpeed(right_, pStar_, rhoStarR_);
            return (S <= shockSpeed) ? starState(right_, rhoStarR_, uStar_) : right_;
        }
        // Right rarefaction.
        const T aR = soundSpeed(right_);
        const T SHead = right_[1] + aR;
        if (S >= SHead)
        {
            return right_;          // untouched right state
        }
        const T STail = uStar_ + soundSpeedStarRight();
        if (S <= STail)
        {
            return starState(right_, rhoStarR_, uStar_);
        }
        const T u = static_cast<T>(2) / (gamma_ + static_cast<T>(1)) *
                    (-aR + static_cast<T>(0.5) * (gamma_ - static_cast<T>(1)) * right_[1] + S);
        const T a = static_cast<T>(2) / (gamma_ + static_cast<T>(1)) *
                    (aR - static_cast<T>(0.5) * (gamma_ - static_cast<T>(1)) * (right_[1] - S));
        const T rho = right_[0] * std::pow(a / aR, static_cast<T>(2) / (gamma_ - static_cast<T>(1)));
        return Vec<T, 3>{rho, u, right_[2] * std::pow(a / aR, static_cast<T>(2) * gamma_ / (gamma_ - static_cast<T>(1)))};
    }

    [[nodiscard]] Vec<T, 3> sampleConserved(const T x, const T t) const
    {
        const Vec<T, 3> W = samplePrimitive(x, t);
        const T rho = W[0];
        const T u = W[1];
        const T p = W[2];
        return Vec<T, 3>{rho, rho * u, p / (gamma_ - static_cast<T>(1)) + static_cast<T>(0.5) * rho * u * u};
    }

    /* Wave speed of the leftmost wave (rarefaction head or left shock). */
    [[nodiscard]] T headSpeed() const
    {
        if (pStar_ > left_[2])
        {
            return shockSpeed(left_, pStar_, rhoStarL_);
        }
        return left_[1] - soundSpeed(left_);
    }

    /* Wave speed of the right shock, or the right rarefaction head. */
    [[nodiscard]] T tailSpeed() const
    {
        if (pStar_ > right_[2])
        {
            return shockSpeed(right_, pStar_, rhoStarR_);
        }
        return right_[1] + soundSpeed(right_);
    }

private:
    T gamma_;
    T x0_;
    Vec<T, 3> left_;
    Vec<T, 3> right_;

    bool converged_ = false;
    T pStar_ = 0;
    T uStar_ = 0;
    T rhoStarL_ = 0;
    T rhoStarR_ = 0;

    [[nodiscard]] T soundSpeed(const Vec<T, 3>& W) const
    {
        return std::sqrt(std::max(gamma_ * W[2] / safeRho(W[0]), static_cast<T>(EPSILON2)));
    }

    /* Shock (Hugoniot) density ratio for a given star pressure. */
    [[nodiscard]] T shockDensity(const Vec<T, 3>& W, const T p) const
    {
        const T ratio = p / std::max(W[2], static_cast<T>(EPSILON2));
        const T numerator = ratio + (gamma_ - static_cast<T>(1)) / (gamma_ + static_cast<T>(1));
        const T denominator = static_cast<T>(1) + (gamma_ - static_cast<T>(1)) / (gamma_ + static_cast<T>(1)) * ratio;
        return W[0] * numerator / denominator;
    }

    /* Star density behind a rarefaction, from the isentropic relation. */
    [[nodiscard]] T rarefactionDensity(const Vec<T, 3>& W, const T p) const
    {
        return W[0] * std::pow(std::max(p / std::max(W[2], static_cast<T>(EPSILON2)), static_cast<T>(EPSILON2)),
                               static_cast<T>(1) / gamma_);
    }

    /* RH shock speed; see developerReference.md before changing. */
    [[nodiscard]] T shockSpeed(const Vec<T, 3>& W, const T p, const T rhoStar) const
    {
        const T deltaP = p - W[2];
        const T rhoK = safeRho(W[0]);
        if (deltaP <= static_cast<T>(0) || rhoStar <= rhoK)
        {
            return W[1]; // not a compression; the caller must use the rarefaction form
        }
        const T denominator = rhoK - rhoK * rhoK / rhoStar;
        if (denominator <= static_cast<T>(EPSILON2))
        {
            return W[1];
        }
        return W[1] + std::sqrt(deltaP / denominator);
    }

    /* Pressure function f_K(p); f_L + f_R = 0 gives the star pressure. */
    void pressureFunction(const Vec<T, 3>& W, const T p, T& f, T& df) const
    {
        const T pSafe = std::max(W[2], static_cast<T>(EPSILON2));
        if (p > pSafe)
        {
            // Shock branch (Rankine-Hugoniot).
            const T coefficient = static_cast<T>(2) / ((gamma_ + static_cast<T>(1)) * safeRho(W[0]));
            const T b = (gamma_ - static_cast<T>(1)) / (gamma_ + static_cast<T>(1));
            const T denominator = p + b * pSafe;
            const T root = std::sqrt(coefficient / denominator);
            f = (p - pSafe) * root;
            df = root * (static_cast<T>(1) - static_cast<T>(0.5) * (p - pSafe) / denominator);
        }
        else
        {
            // Rarefaction branch (isentropic).
            const T aK = soundSpeed(W);
            const T pRatio = std::max(p / pSafe, static_cast<T>(EPSILON2));
            const T power = std::pow(pRatio, (gamma_ - static_cast<T>(1)) / (static_cast<T>(2) * gamma_));
            f = static_cast<T>(2) * aK / (gamma_ - static_cast<T>(1)) * (power - static_cast<T>(1));
            df = power / (aK * p);
        }
    }

    void solve()
    {
        const T aL = soundSpeed(left_);
        const T aR = soundSpeed(right_);

        // Initial pressure guess. Both branch markers below only choose a seed for the
        // root find; correctness does not depend on the seed.
        const T massFluxGuess = static_cast<T>(0.5) * (left_[0] + right_[0]) * static_cast<T>(0.5) * (aL + aR);
        const T pPV = std::max(static_cast<T>(0),
                               static_cast<T>(0.5) * (left_[2] + right_[2]) -
                               static_cast<T>(0.5) * (right_[1] - left_[1]) * massFluxGuess);
        const T pMin = std::min(left_[2], right_[2]);
        const T pMax = std::max(left_[2], right_[2]);
        T p;
        if (pPV >= pMin)
        {
            // Two-rarefaction guess.
            const T gamma1 = static_cast<T>(0.5) * (gamma_ - static_cast<T>(1)) / gamma_;
            p = std::pow((aL + aR - static_cast<T>(0.5) * (gamma_ - static_cast<T>(1)) * (right_[1] - left_[1])) /
                         (aL / std::pow(left_[2], gamma1) + aR / std::pow(right_[2], gamma1)),
                         static_cast<T>(1) / gamma1);
            p = std::max(p, static_cast<T>(EPSILON2));
        }
        else
        {
            // Two-shock guess.
            const T gL = std::sqrt(static_cast<T>(2) / ((gamma_ + static_cast<T>(1)) * safeRho(left_[0])));
            const T gR = std::sqrt(static_cast<T>(2) / ((gamma_ + static_cast<T>(1)) * safeRho(right_[0])));
            p = std::max(static_cast<T>(EPSILON2),
                         (gL * left_[2] + gR * right_[2] - (right_[1] - left_[1])) / (gL + gR));
        }

        T fL = 0, fR = 0, dfL = 0, dfR = 0;
        auto F = [&](const T pressure)
        {
            pressureFunction(left_, pressure, fL, dfL);
            pressureFunction(right_, pressure, fR, dfR);
            return fL + fR;
        };

        // Establish a bracket. The seeds cover both extremes (strong rarefaction,
        // strong shock) and are widened geometrically until F changes sign.
        T pLower = static_cast<T>(1e-10);
        T pUpper = std::max(pMax, p);
        T fLower = F(pLower);
        T fUpper = F(pUpper);
        for (uint widen = 0; widen < 200 && (fLower > static_cast<T>(0) || fUpper < static_cast<T>(0)); ++widen)
        {
            if (fLower > static_cast<T>(0)) pLower *= static_cast<T>(0.1);
            if (fUpper < static_cast<T>(0)) pUpper *= static_cast<T>(10);
            pLower = std::max(pLower, static_cast<T>(1e-30));
            pUpper = std::min(pUpper, static_cast<T>(1e30));
            fLower = F(pLower);
            fUpper = F(pUpper);
        }

        p = std::clamp(p, pLower, pUpper);
        const T tolerance = static_cast<T>(1e-14) * std::max(static_cast<T>(1), pMax);
        bool bracketed = (fLower <= static_cast<T>(0) && fUpper >= static_cast<T>(0));
        if (bracketed)
        {
            for (uint iteration = 0; iteration < 200; ++iteration)
            {
                const T fMid = F(p);
                if (std::abs(fMid) <= tolerance || (pUpper - pLower) <= tolerance * static_cast<T>(1e-2))
                {
                    break;
                }
                if (fMid > static_cast<T>(0))
                {
                    pUpper = p;
                    fUpper = fMid;
                    // Illinois: halve the retained end's magnitude to keep superlinear
                    // convergence without losing the bracket.
                    fLower *= static_cast<T>(0.5);
                }
                else
                {
                    pLower = p;
                    fLower = fMid;
                    fUpper *= static_cast<T>(0.5);
                }
                // Regula falsi step on the (possibly modified) bracket.
                const T denominator = fUpper - fLower;
                const T pFalse = (std::abs(denominator) > static_cast<T>(EPSILON3))
                                     ? (pLower * fUpper - pUpper * fLower) / denominator
                                     : static_cast<T>(0.5) * (pLower + pUpper);
                p = (pFalse > pLower && pFalse < pUpper) ? pFalse : static_cast<T>(0.5) * (pLower + pUpper);
            }
        }
        else
        {
            // Should be unreachable for physically realisable data; fall back to the
            // midpoint so the result is at least finite and reported as unconverged.
            p = static_cast<T>(0.5) * (pLower + pUpper);
        }

        pStar_ = p;
        pressureFunction(left_, pStar_, fL, dfL);
        pressureFunction(right_, pStar_, fR, dfR);
        const T uFromLeft = left_[1] - fL;
        const T uFromRight = right_[1] + fR;
        uStar_ = static_cast<T>(0.5) * (uFromLeft + uFromRight);

        // Star densities behind the left/right waves.
        rhoStarL_ = (pStar_ > left_[2]) ? shockDensity(left_, pStar_) : rarefactionDensity(left_, pStar_);
        rhoStarR_ = (pStar_ > right_[2]) ? shockDensity(right_, pStar_) : rarefactionDensity(right_, pStar_);

        // A converged solution has consistent star velocities from both sides and a
        // positive post-wave density.
        converged_ = std::isfinite(pStar_) && std::isfinite(uStar_) &&
                     pStar_ > static_cast<T>(0) &&
                     rhoStarL_ > static_cast<T>(0) && rhoStarR_ > static_cast<T>(0) &&
                     std::abs(uFromLeft - uFromRight) <=
                         static_cast<T>(1e-6) * std::max(static_cast<T>(1), std::abs(uStar_)) + static_cast<T>(1e-10);
    }

    [[nodiscard]] Vec<T, 3> starState(const Vec<T, 3>& W, const T rhoStar, const T uStar) const
    {
        (void)W;
        // Pressure is uniform across the contact, so the star pressure applies.
        return Vec<T, 3>{rhoStar, uStar, pStar_};
    }
};


template <typename T>
class SolutionValidator
{
public:
    SolutionValidator(const ExactRiemannSolver<T>& exact, const VariableFormulation formulation)
        : exact_(exact), formulation_(formulation)
    {
    }


    VolumetricData<T> validate(const Mesh<T, 1>& mesh, const T time,
                               const std::vector<FieldSpec>& requestedFields) const
    {
        VolumetricData<T> out;
        out.time = time;
        out.numCells = mesh.numCells;
        out.formulation = formulation_;
        out.exactInfo.valid = exact_.converged();
        out.exactInfo.pStar = exact_.starPressure();
        out.exactInfo.uStar = exact_.starVelocity();
        out.exactInfo.rhoStarL = exact_.starDensityLeft();
        out.exactInfo.rhoStarR = exact_.starDensityRight();
        out.exactInfo.contactSpeed = exact_.starVelocity();
        out.exactInfo.shockSpeed = exact_.tailSpeed();
        out.exactInfo.headSpeed = exact_.headSpeed();
        out.exactInfo.tailSpeed = exact_.tailSpeed();

        const uint n = mesh.numCells;
        if (n == 0)
        {
            return out;
        }

        out.fields.resize(requestedFields.size());
        for (size_t f = 0; f < requestedFields.size(); ++f)
        {
            FieldStats<T>& stats = out.fields[f];
            stats.name = requestedFields[f].name;
            stats.unit = requestedFields[f].unit;
            stats.conserved = requestedFields[f].conserved;
            stats.x.reserve(n);
            stats.exact.reserve(n);
            stats.simulated.reserve(n);
            stats.error.reserve(n);
        }

        std::vector<T> sumAbs(requestedFields.size(), static_cast<T>(0));
        std::vector<T> sumSq(requestedFields.size(), static_cast<T>(0));
        std::vector<T> sumAbsExact(requestedFields.size(), static_cast<T>(0));
        std::vector<T> sumSqExact(requestedFields.size(), static_cast<T>(0));
        std::vector<T> maxAbs(requestedFields.size(), static_cast<T>(0));
        std::vector<T> maxAbsExact(requestedFields.size(), static_cast<T>(0));

        ConservationStats<T>& cons = out.conservation;

        for (uint cell = 0; cell < n; ++cell)
        {
            const T x = mesh.cells.cellCentersPositions[cell][0];

            // Normalise the simulated state to both formulations, so the reported
            // fields are the same regardless of how the caller stores its data.
            const Vec<T, 3>& raw = mesh.cells.cellAverages[cell];
            const Vec<T, 3> WSim = asPrimitive<T, 1>(raw, formulation_);
            const Vec<T, 3> USim = asConserved<T, 1>(raw, formulation_);

            // Exact state at the same location, in the same normalisations.
            const Vec<T, 3> WExact = exact_.samplePrimitive(x, time);
            const Vec<T, 3> UExact = exact_.sampleConserved(x, time);

            for (size_t f = 0; f < requestedFields.size(); ++f)
            {
                const T sim = extractField<T, 1>(requestedFields[f].variable, WSim, USim);
                const T exa = extractField<T, 1>(requestedFields[f].variable, WExact, UExact);
                const T err = sim - exa;

                FieldStats<T>& stats = out.fields[f];
                stats.x.push_back(x);
                stats.simulated.push_back(sim);
                stats.exact.push_back(exa);
                stats.error.push_back(err);

                sumAbs[f] += std::abs(err);
                sumSq[f] += err * err;
                sumAbsExact[f] += std::abs(exa);
                sumSqExact[f] += exa * exa;
                maxAbs[f] = std::max(maxAbs[f], std::abs(err));
                maxAbsExact[f] = std::max(maxAbsExact[f], std::abs(exa));
            }

            cons.massSimulated += USim[0] * mesh.cells.cellVolumes[cell];
            cons.momentumSimulated += USim[1] * mesh.cells.cellVolumes[cell];
            cons.energySimulated += USim[2] * mesh.cells.cellVolumes[cell];

            out.cellCenters.push_back(x);
        }

        // Quadrature intervals as a multiple of 2*(n-1) to land on cell faces.
        const uint baseIntervals = std::max<uint>(64u, 4u * (n > 1 ? n - 1 : 1u));
        const uint numIntervals = 2u * ((baseIntervals + 1u) / 2u);
        const T xLeft = (n > 0) ? mesh.cells.cellCentersPositions[0][0] - static_cast<T>(0.5) * mesh.cells.cellVolumes[0]
                                : mesh.boundingBox.first[0];
        const T xRight = (n > 0) ? mesh.cells.cellCentersPositions[n - 1][0] + static_cast<T>(0.5) * mesh.cells.cellVolumes[n - 1]
                                 : mesh.boundingBox.second[0];


        const T interfacePosition = exact_.interfacePosition();
        const T leftSpan = std::clamp(interfacePosition, xLeft, xRight) - xLeft;
        const T rightSpan = xRight - std::clamp(interfacePosition, xLeft, xRight);
        const Vec<T, 3> initialLeft = asConserved<T, 1>(exact_.leftState(), VariableFormulation::PRIMITIVE);
        const Vec<T, 3> initialRight = asConserved<T, 1>(exact_.rightState(), VariableFormulation::PRIMITIVE);
        cons.massExact = initialLeft[0] * leftSpan + initialRight[0] * rightSpan;
        cons.momentumExact = initialLeft[1] * leftSpan + initialRight[1] * rightSpan;
        cons.energyExact = initialLeft[2] * leftSpan + initialRight[2] * rightSpan;


        for (size_t f = 0; f < requestedFields.size(); ++f)
        {
            const FieldVariable variable = requestedFields[f].variable;
            if (variable == FieldVariable::DENSITY)
            {
                out.fields[f].exactIntegral = cons.massExact;
            }
            else if (variable == FieldVariable::MOMENTUM)
            {
                out.fields[f].exactIntegral = cons.momentumExact;
            }
            else if (variable == FieldVariable::ENERGY)
            {
                out.fields[f].exactIntegral = cons.energyExact;
            }
            else
            {
                out.fields[f].exactIntegral = integrateExactField(variable, xLeft, xRight, numIntervals, time);
            }
        }

        // Per-field pointwise statistics and the signed integral error against the
        // reference computed above.
        for (size_t f = 0; f < requestedFields.size(); ++f)
        {
            FieldStats<T>& stats = out.fields[f];
            const T invN = static_cast<T>(1) / static_cast<T>(n);
            stats.l1 = sumAbs[f] * invN;
            stats.l2 = std::sqrt(sumSq[f] * invN);
            stats.lInf = maxAbs[f];
            const T meanAbsExact = sumAbsExact[f] * invN;
            const T meanSqExact = sumSqExact[f] * invN;
            stats.relL1 = (meanAbsExact > static_cast<T>(EPSILON2)) ? stats.l1 / meanAbsExact : static_cast<T>(0);
            stats.relL2 = (meanSqExact > static_cast<T>(EPSILON2)) ? stats.l2 / std::sqrt(meanSqExact) : static_cast<T>(0);
            stats.relLInf = (maxAbsExact[f] > static_cast<T>(EPSILON2)) ? stats.lInf / maxAbsExact[f] : static_cast<T>(0);

            T simulatedIntegral = static_cast<T>(0);
            for (uint cell = 0; cell < n; ++cell)
            {
                simulatedIntegral += stats.simulated[cell] * mesh.cells.cellVolumes[cell];
            }
            stats.errorIntegral = simulatedIntegral - stats.exactIntegral;
        }

        auto relative = [](const T error, const T reference)
        {
            return (std::abs(reference) > static_cast<T>(EPSILON2)) ? error / reference : static_cast<T>(0);
        };
        cons.massError = cons.massSimulated - cons.massExact;
        cons.momentumError = cons.momentumSimulated - cons.momentumExact;
        cons.energyError = cons.energySimulated - cons.energyExact;
        cons.massErrorRel = relative(cons.massError, cons.massExact);
        cons.momentumErrorRel = relative(cons.momentumError, cons.momentumExact);
        cons.energyErrorRel = relative(cons.energyError, cons.energyExact);

        return out;
    }

    /* Compare with the default field set for this formulation. */
    VolumetricData<T> validate(const Mesh<T, 1>& mesh, const T time) const
    {
        return validate(mesh, time, defaultFields(formulation_));
    }


    static void printReport(const VolumetricData<T>& data)
    {
        std::printf("Exact-solution comparison\n");
        std::printf("  formulation : %s\n",
                    data.formulation == VariableFormulation::CONSERVATIVE ? "conservative" : "primitive");
        std::printf("  time        : %.6g\n", static_cast<double>(data.time));
        std::printf("  cells       : %u\n", data.numCells);
        if (data.exactInfo.valid)
        {
            std::printf("  exact p*    : %.6g   u* : %.6g\n",
                        static_cast<double>(data.exactInfo.pStar), static_cast<double>(data.exactInfo.uStar));
            std::printf("  shock speed : %.6g   contact : %.6g\n",
                        static_cast<double>(data.exactInfo.shockSpeed), static_cast<double>(data.exactInfo.contactSpeed));
        }
        else
        {
            std::printf("  exact solver did NOT converge; values below are unreliable\n");
        }
        std::printf("  %-8s %12s %12s %12s %12s %12s\n", "field", "L1", "L2", "Linf", "relL1", "int.error");
        for (const FieldStats<T>& field : data.fields)
        {
            std::printf("  %-8s %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                        field.name.c_str(),
                        static_cast<double>(field.l1), static_cast<double>(field.l2),
                        static_cast<double>(field.lInf), static_cast<double>(field.relL1),
                        static_cast<double>(field.errorIntegral));
        }
        std::printf("  conservation: mass %.3e (rel %.3e), momentum %.3e (rel %.3e), energy %.3e (rel %.3e)\n",
                    static_cast<double>(data.conservation.massError), static_cast<double>(data.conservation.massErrorRel),
                    static_cast<double>(data.conservation.momentumError), static_cast<double>(data.conservation.momentumErrorRel),
                    static_cast<double>(data.conservation.energyError), static_cast<double>(data.conservation.energyErrorRel));
    }

private:
    /* Composite Simpson rule on a uniform grid. */
    static T simpson(const std::function<T(T)>& f, const T a, const T b, uint numIntervals)
    {
        if (numIntervals % 2 != 0) ++numIntervals;
        if (numIntervals == 0 || b <= a)
        {
            return static_cast<T>(0);
        }
        const T h = (b - a) / static_cast<T>(numIntervals);
        T total = f(a) + f(b);
        for (uint i = 1; i < numIntervals; ++i)
        {
            total += ((i % 2 == 0) ? static_cast<T>(2) : static_cast<T>(4)) * f(a + static_cast<T>(i) * h);
        }
        return total * h / static_cast<T>(3);
    }


    T integrateExactField(const FieldVariable field, const T a, const T b, const uint numIntervals, const T time) const
    {
        if (b <= a || time <= static_cast<T>(0) || !exact_.converged())
        {
            return static_cast<T>(0);
        }
        auto sample = [this, field, time](const T x)
        {
            const Vec<T, 3> W = exact_.samplePrimitive(x, time);
            return extractField<T, 1>(field, W, exact_.sampleConserved(x, time));
        };

        const T x0 = exact_.interfacePosition();
        const T xHead = x0 + time * exact_.headSpeed();
        const T xShock = x0 + time * exact_.tailSpeed();
        const T xContact = x0 + time * exact_.starVelocity();
        // The rarefaction tail only exists when the left wave is a rarefaction.
        const bool leftRarefaction = exact_.starPressure() < exact_.leftState()[2];
        const T xTail = leftRarefaction ? x0 + time * (exact_.starVelocity() - exact_.soundSpeedStarLeft())
                                        : xContact;

        std::vector<T> cuts;
        cuts.push_back(a);
        for (const T candidate : {xHead, xTail, xContact, xShock})
        {
            if (candidate > a && candidate < b)
            {
                cuts.push_back(candidate);
            }
        }
        cuts.push_back(b);
        std::sort(cuts.begin(), cuts.end());

        T total = 0;
        for (size_t i = 1; i < cuts.size(); ++i)
        {
            if (cuts[i] - cuts[i - 1] <= static_cast<T>(0))
            {
                continue;
            }
            // Enough sub-intervals per region to resolve the fan.
            // The constant regions are integrated exactly by any rule.
            total += simpson(sample, cuts[i - 1], cuts[i], std::max<uint>(numIntervals, 64u));
        }
        return total;
    }

    const ExactRiemannSolver<T>& exact_;
    VariableFormulation formulation_;
};

} // namespace riemann
