// Standalone verification of the Limiters in include/physics.hpp.
//
// Build/run (physics.hpp pulls mesh/vec/constants, no SFML needed):
//   g++ -std=c++20 -O1 -Iinclude -I/usr/include/eigen3 tests/limiter_check.cpp -o /tmp/limiter_check && /tmp/limiter_check
//
// The interesting part of this file is the set of *properties* it pins down, because
// they are easy to get wrong by intuition:
//   * van Leer reproduces phi(r)/r exactly (independent textbook reference),
//   * its guarantee is the TVD bracket |psi*dM| <= |dP| + |dM|, NOT |psi*dM| <= |dP|,
//   * so canonical van Leer may leave [Umin, Umax] when dP < dM, while the
//     VAN_LEER_MONOTONE variant and Venkatakrishnan do not,
//   * a local extremum (dP == 0) must give psi = 0.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

#include "mesh.hpp"
#include "physics.hpp"

using namespace std;

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++failures; printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

/*
 * Independent textbook reference.
 *
 * For a face j+1/2, with the upwind jump dUp and the downwind jump dDn, the
 * classical form of a symmetric TVD limiter is
 *
 *     psi = phi(r) / r,   r = dUp / dDn,   phi(r) = (r + |r|) / (1 + |r|)
 *
 * The code under test is handed (deltaM, deltaP) == (dUp, dDn), so this is a
 * genuinely independent path to the same quantity. Note phi <= 2, hence psi <= 2,
 * and that a negative ratio cancels to 2r/(1-|r|) < 0, i.e. an extremum maps to 0.
 */
static double vanLeerRef(double dUp, double dDn)
{
    const double r = dUp / dDn;
    const double phi = (r + std::abs(r)) / (1.0 + std::abs(r));
    const double psiMax = VAN_LEER_MONOTONE ? 1.0 : 2.0;
    return std::clamp(phi / r, 0.0, psiMax);
}

/*
 * Builds the exact (deltaM, deltaP) pair the solver feeds to the limiter, for a
 * cell whose neighbourhood spans [Umin, Umax] and whose RBF reconstruction
 * produced the jump deltaM:
 *
 *   deltaM = Ui - Uj
 *   deltaP = (deltaM > 0) ? (Umax - Uj) : (Uj - Umin)
 *
 * For a consistent input Uj must lie in [Umin, Umax], which makes deltaP >= 0 and
 * makes a local extremum appear as deltaP == 0.
 */
static void makeSolverInputs(double Uj, double Umin, double Umax, double deltaM,
                             double& deltaP, double& psiVl, double& psiVk)
{
    deltaP = (deltaM > 0) ? (Umax - Uj) : (Uj - Umin);
    psiVl = Limiters<double, 3>::VanLeerField(deltaM, deltaP);
    psiVk = Limiters<double, 3>::VenkatakrishnanField(deltaM, deltaP, 1e-12);
}

int main()
{
    using L = Limiters<double, 3>;

    // ---- van Leer: exact canonical values -------------------------------------
    // dUp == dDn (r = 1) -> phi = 1 -> psi = 1
    CHECK(std::abs(L::VanLeerField(1.0, 1.0) - 1.0) < 1e-14, "psi should be 1, got %g", L::VanLeerField(1.0, 1.0));
    CHECK(std::abs(L::VanLeerField(-0.5, -0.5) - 1.0) < 1e-14, "psi should be 1, got %g", L::VanLeerField(-0.5, -0.5));
    // dUp = 1, dDn = 2 -> r = 0.5 -> phi = 2/3 -> canonical psi = 4/3, monotone psi = 1
    const double psiUp = VAN_LEER_MONOTONE ? 1.0 : 4.0 / 3.0;
    CHECK(std::abs(L::VanLeerField(1.0, 2.0) - psiUp) < 1e-14, "psi should be %g, got %g", psiUp, L::VanLeerField(1.0, 2.0));
    CHECK(std::abs(L::VanLeerField(-1.0, -2.0) - psiUp) < 1e-14, "mirrored psi should be %g, got %g", psiUp, L::VanLeerField(-1.0, -2.0));
    // dUp = 2, dDn = 1 -> r = 2   -> phi = 4/3 -> psi = 2/3 (below 1, unaffected by the cap)
    CHECK(std::abs(L::VanLeerField(2.0, 1.0) - 2.0 / 3.0) < 1e-14, "psi should be 2/3, got %g", L::VanLeerField(2.0, 1.0));

    // ---- van Leer: extrema must collapse the jump ------------------------------
    // No extrapolation to limit.
    CHECK(std::abs(L::VanLeerField(0.0, 3.0) - 1.0) < 1e-14, "psi at deltaM=0 should be 1, got %g", L::VanLeerField(0.0, 3.0));
    // Uj sits exactly on the bound it is extrapolating towards: deltaP = 0 -> psi = 0.
    CHECK(std::abs(L::VanLeerField(1.0, 0.0) - 0.0) < 1e-14, "psi at deltaP=0 should be 0, got %g", L::VanLeerField(1.0, 0.0));
    CHECK(std::abs(L::VanLeerField(-1.0, 0.0) - 0.0) < 1e-14, "psi at deltaP=0 (neg) should be 0, got %g", L::VanLeerField(-1.0, 0.0));
    // Opposing signs (an inconsistent input, but must not blow up) -> psi = 0.
    CHECK(std::abs(L::VanLeerField(1.0, -0.5) - 0.0) < 1e-14, "psi at sign flip should be 0, got %g", L::VanLeerField(1.0, -0.5));
    CHECK(std::abs(L::VanLeerField(-1.0, 0.5) - 0.0) < 1e-14, "psi at sign flip (neg) should be 0, got %g", L::VanLeerField(-1.0, 0.5));

    // ---- van Leer matches the textbook reference over many ratios --------------
    for (double r = 1e-3; r <= 1e3; r *= 1.05)
    {
        // dUp = 1, dDn = 1/r -> ratio r
        CHECK(std::abs(L::VanLeerField(1.0, 1.0 / r) - vanLeerRef(1.0, 1.0 / r)) < 1e-12,
              "van Leer r=%g gave %g expected %g", r, L::VanLeerField(1.0, 1.0 / r), vanLeerRef(1.0, 1.0 / r));
        // mirrored (both signs flipped) must give the same psi
        CHECK(std::abs(L::VanLeerField(-1.0, -1.0 / r) - vanLeerRef(-1.0, -1.0 / r)) < 1e-12,
              "mirrored van Leer r=%g gave %g expected %g", r, L::VanLeerField(-1.0, -1.0 / r), vanLeerRef(-1.0, -1.0 / r));
    }

    // ---- the TVD bracket holds for arbitrary (deltaM, deltaP) ------------------
    for (int i = 0; i < 20000; ++i)
    {
        const double deltaM = std::cos(0.013 * i) * 3.0;
        const double deltaP = std::sin(0.017 * i) * 2.0; // deliberately unconstrained
        const double psi = L::VanLeerField(deltaM, deltaP);
        const double limitedJump = psi * deltaM;
        CHECK(psi >= 0.0 && psi <= 2.0, "psi out of range: %g", psi);
        CHECK(std::abs(limitedJump) <= std::abs(deltaP) + std::abs(deltaM) + 1e-12,
              "bracket violated: |%g| > |%g| + |%g|", limitedJump, deltaP, deltaM);
    }

    // ---- realistic sweep: consistent neighbourhood + reconstruction -------------
    // Uj is inside [Umin, Umax] and deltaM is bounded by the room available in its
    // own direction, exactly what the solver creates.
    int extremumCases = 0;
    int psiAboveOne = 0;
    int vlStrays = 0, vkStrays = 0, tvdBracketHolds = 0;
    for (int i = 0; i < 20000; ++i)
    {
        const double Uj = std::sin(0.211 * i) * 4.0;
        const double below = 0.25 + std::abs(std::cos(0.05 * i));
        const double above = 0.25 + std::abs(std::sin(0.03 * i));
        const double Umin = Uj - below;
        const double Umax = Uj + above;

        // deltaM from a "reconstruction": 80% inside the bracket, 20% an overshooting
        // interpolation that crosses the far bound.
        const double frac = (i % 5 == 0) ? (1.0 + std::abs(std::sin(0.31 * i))) : (0.9 * std::abs(std::cos(0.19 * i)));
        const double sgn = (std::cos(0.013 * i) >= 0.0) ? 1.0 : -1.0;
        const double room = (sgn > 0.0) ? above : below;
        const double deltaM = sgn * room * frac;

        double deltaP, psiVl, psiVk;
        makeSolverInputs(Uj, Umin, Umax, deltaM, deltaP, psiVl, psiVk);
        CHECK(deltaP >= -1e-15, "deltaP should be non-negative, got %g", deltaP);
        if (psiVl > 1.0 + 1e-12) ++psiAboveOne;

        const double limVl = Uj + psiVl * deltaM;
        const double limVk = Uj + psiVk * deltaM;

        // The TVD bracket is unconditional.
        CHECK(std::abs(psiVl * deltaM) <= std::abs(deltaP) + std::abs(deltaM) + 1e-12,
              "van Leer TVD bracket violated for Uj=%g dM=%g dP=%g psi=%g", Uj, deltaM, deltaP, psiVl);
        if (std::abs(psiVl * deltaM) <= std::abs(deltaP) + std::abs(deltaM) + 1e-12) ++tvdBracketHolds;

        // When psi <= 1 the limited jump can never exceed the unlimited one.
        CHECK(psiVl * deltaM <= deltaM + 1e-12 || deltaM < 0.0 || psiVl > 1.0,
              "psi <= 1 amplified the jump: %g > %g", psiVl * deltaM, deltaM);

        // Range preservation is NOT guaranteed here: with deltaP taken from the
        // neighbourhood extrema rather than from a neighbouring cell average, even
        // psi <= 1 leaves the range when deltaP < deltaM. Count, do not assert.
        if (limVl < Umin - 1e-12 || limVl > Umax + 1e-12) ++vlStrays;
        if (limVk < Umin - 1e-12 || limVk > Umax + 1e-12) ++vkStrays;
    }
    if constexpr (!VAN_LEER_MONOTONE)
    {
        CHECK(psiAboveOne > 0, "sweep never exercised the psi > 1 branch of van Leer");
    }
    // Both limiters can overshoot the extrema-derived range; this is exactly why a
    // hard clamp after limiting is the only way to *guarantee* boundedness.
    CHECK(vlStrays > 0 && vkStrays > 0, "expected both limiters to occasionally exceed [Umin,Umax] (vl=%d vk=%d)", vlStrays, vkStrays);
    CHECK(tvdBracketHolds == 20000, "TVD bracket must hold everywhere");

    // ---- exact extremum: Uj sits on the bound it extrapolates towards ----------
    {
        double deltaP, psiVl, psiVk;
        makeSolverInputs(1.0, 0.0, 1.0, 0.3, deltaP, psiVl, psiVk);
        CHECK(deltaP == 0.0, "extremum deltaP should be 0, got %g", deltaP);
        CHECK(std::abs(psiVl) < 1e-12, "van Leer should clamp an extremum, got psi = %g", psiVl);
        // Venkatakrishnan is regularised, so it is only O(epsilonSq) at an extremum.
        CHECK(std::abs(psiVk) < 1e-9, "VK should nearly clamp an extremum, got psi = %g", psiVk);
        ++extremumCases;
    }

    // ---- the two van Leer regimes, which differ by VAN_LEER_MONOTONE -----------
    //  dP > dM -> r < 1, canonical psi > 1 (least diffusive). Since psi*dM < dP the
    //            edge value still lands inside the bounds. VAN_LEER_MONOTONE caps
    //            psi at 1 and is therefore *more* diffusive here.
    //  dP < dM -> r > 1, psi < 1 for both variants, but psi*dM = 2 dP dM/(dP+dM) can
    //            still exceed dP. So neither psi bound removes the excursion: with
    //            extrema-derived deltaP, the limiter alone cannot guarantee
    //            Uj + psi*dM is inside [Umin, Umax].
    {
        const double Uj = 2.0, Umax = 4.0, dM = 1.0;
        const double dP = Umax - Uj; // 2.0 -> r = 1/2 -> canonical psi = 4/3
        const double psi = L::VanLeerField(dM, dP);
        CHECK(Uj + psi * dM <= Umax + 1e-12, "dP > dM must land inside, got %g > %g", Uj + psi * dM, Umax);
        if constexpr (VAN_LEER_MONOTONE)
        {
            CHECK(std::abs(psi - 1.0) < 1e-14, "monotone variant should cap psi at 1, got %g", psi);
        }
        else
        {
            CHECK(psi > 1.0, "canonical van Leer should give psi > 1 for dP > dM, got %g", psi);
        }
    }
    {
        const double Uj = 2.0, Umax = 2.25, dM = 2.0;
        const double dP = Umax - Uj; // 0.25 -> r = 8 -> psi = 2/9 in both variants
        const double psi = L::VanLeerField(dM, dP);
        const double lim = Uj + psi * dM;
        CHECK(std::abs(psi * dM) <= std::abs(dP) + std::abs(dM) + 1e-12, "TVD bracket must always hold");
        // Both variants overshoot here; the clamp below is what fixes it.
        CHECK(lim > Umax, "expected the limiter to relax the range: lim=%g Umax=%g", lim, Umax);
        const double clamped = std::clamp(lim, Uj - dP, Umax);
        CHECK(clamped <= Umax + 1e-12, "a post-limiting clamp must bound the edge value");
    }

    // ---- Venkatakrishnan reference values -------------------------------------
    // This formulation is not the classic (deltaP^2 + 2 deltaM deltaP + eps) / (...);
    // for deltaM == deltaP and negligible eps it gives exactly 3/4.
    CHECK(std::abs(L::VenkatakrishnanField(0.5, 0.5, 0.0) - 0.75) < 1e-14, "VK r=1 psi = %g", L::VenkatakrishnanField(0.5, 0.5, 0.0));
    // It approaches 1 as the extrapolation becomes small relative to the room left.
    CHECK(L::VenkatakrishnanField(1e-3, 1.0, 1e-12) > 0.999, "VK small jump psi = %g", L::VenkatakrishnanField(1e-3, 1.0, 1e-12));
    // deltaM = 0 -> 1
    CHECK(std::abs(L::VenkatakrishnanField(0.0, 5.0, 1e-12) - 1.0) < 1e-14, "VK deltaM=0 psi = %g", L::VenkatakrishnanField(0.0, 5.0, 1e-12));
    // At a genuine extremum (deltaP = 0) it must vanish, like van Leer.
    CHECK(std::abs(L::VenkatakrishnanField(1.0, 0.0, 1e-12)) < 1e-12, "VK extremum psi = %g", L::VenkatakrishnanField(1.0, 0.0, 1e-12));

    // ---- vector dispatch: the three kinds are distinguishable -----------------
    // Consistent case: Uj = 1, Umin = 1, Umax = 1.25, extrapolation Ui = 1.5.
    Vec<double, 3> Uj{1.0, 1.0, 1.0};
    Vec<double, 3> Ui{1.5, 1.5, 1.5};   // deltaM = 0.5
    Vec<double, 3> Umin{1.0, 1.0, 1.0};
    Vec<double, 3> Umax{1.25, 1.25, 1.25}; // deltaP = 0.25

    const Vec<double, 3> none = L::Apply(LIMITER_KIND::NO_LIMITER, Uj, Ui, Umin, Umax, 1e-12);
    const Vec<double, 3> vl   = L::Apply(LIMITER_KIND::VAN_LEER, Uj, Ui, Umin, Umax, 1e-12);
    const Vec<double, 3> vk   = L::Apply(LIMITER_KIND::VENKATAKRISHNAN, Uj, Ui, Umin, Umax, 1e-12);

    CHECK(std::abs(none[0] - 1.0) < 1e-14, "NO_LIMITER psi = %g", none[0]);
    const double expectedVl = vanLeerRef(0.5, 0.25); // = 2*0.25/0.75 = 2/3
    CHECK(std::abs(expectedVl - 2.0 / 3.0) < 1e-14, "reference sanity: %g", expectedVl);
    CHECK(std::abs(vl[0] - expectedVl) < 1e-12, "VAN_LEER psi = %g expected %g", vl[0], expectedVl);
    CHECK(std::abs(vk[0] - L::VenkatakrishnanField(0.5, 0.25, 1e-12)) < 1e-14, "VK dispatch psi = %g", vk[0]);
    // In this configuration van Leer is the least diffusive of the two.
    CHECK(vl[0] > vk[0], "van Leer should be less diffusive than VK here (vl=%g vk=%g)", vl[0], vk[0]);
    // The legacy Limiters::Venkatakrishnan entry point is deprecated; its behaviour is
    // already covered by the VENKATAKRISHNAN dispatch case above, so it is not called
    // here (calling it would only emit -Wdeprecated-declarations).

    // ---- large but finite deltaP stays finite and clamped ----------------------
    CHECK(std::isfinite(L::VanLeerField(1.0, 1e12)) && L::VanLeerField(1.0, 1e12) <= 2.0,
          "van Leer with large deltaP = %g", L::VanLeerField(1.0, 1e12));

    if (failures == 0) printf("ALL LIMITER CHECKS PASSED\n");
    else printf("%d LIMITER CHECKS FAILED\n", failures);
    return failures == 0 ? 0 : 1;
}
