
TODO:

1. Integration
Numerical integration is in place (composite 4 x 16 Gauss-Legendre). For 2D, only the
quadrature domain changes.

We first get the rbfs-cell averages relations of the form:
u = Ax, where u are the vector of cell averages, x is the value of the rbs, A is the coefficients of interactions thata ccounts for the positions of the RBFs in the domain.

A is constant throughout the simulation since the RBFs' centers do not move. it would thus be worthwhile to take an inverse of the matrix rather than use an iterative method. 

Alternatively, we could find partial(x)/partial(u). Since it is linear, I just need to do 1 base computation to find x, then do N calculations of u_modified, where N is the number of cells (also the dimension of the x and u vector). u_modified (s) will be computed via just incrementing/decrementing u in one place by a small amount.






## Compact RBF reconstruction

Kernels: Wendland C0/C2/C4/C6 and D4C4, Wu C2/C5, compact exponential/rational, and a
C-infinity bump. Limiters: Venkatakrishnan, van Leer, minmod, superbee, MC, van Albada,
Barth-Jespersen, Koren. The app has kernel, support, limiter and bracket controls.

The local cell-average RBF fit now carries constant and linear polynomial
reproduction, plus a 1D quadratic term, and the edge reconstruction is a
mean-preserving quadratic that is third order on smooth data. The limiter is applied
to primitive `{rho, u, p}` extrema, which keeps the star-region velocity/pressure
oscillation under control. See [`rbfThirdOrderStudy.md`](rbfThirdOrderStudy.md),
[`rbfBoundaryDipStudy.md`](rbfBoundaryDipStudy.md) and
[`rbfX5UndershootStudy.md`](rbfX5UndershootStudy.md).

```bash
cmake -S . -B cmake-build-debug -DRBF_BUILD_TOOLS=ON
cmake --build cmake-build-debug --target kernel_check limiter_check rbf_kernel_study -j"$(nproc)"
./cmake-build-debug/kernel_check && ./cmake-build-debug/limiter_check
./cmake-build-debug/rbf_kernel_study --sod-cells 100,200 --rec-cells 50,100 \
    --tend 2.0 --cfl 0.0125 --support-sweep 0.4,1.9 --kernels all \
    --limiters venkatakrishnan --sod-only --out results/kernel_classes
```

Reports: [`rbfKernelStudy.md`](rbfKernelStudy.md),
[`rbfKernelLimiterStudy.md`](rbfKernelLimiterStudy.md). Raw results are under
`results/` (gitignored).

## Sod benchmark validation

`data/sod_shock_tube.csv` is the NASA NPARC Alliance shock-tube reference profile
(`stp.dx`, `str.dx`, `stu.dx`, `stm.dx` from the
[NPARC validation archive](https://www.grc.nasa.gov/WWW/wind/valid/stube/stube.html)).
`include/sod_shock_tube.hpp` loads the CSV and its metadata, `include/sod_simulation.hpp`
runs the 1D Euler solver with the standard Sod states and transmissive ends, and
`include/sod_validation.hpp` compares the cell solution at the published data points.
The app no longer uses the exact Riemann solver for validation.

```bash
cmake -S . -B cmake-build-debug -DRBF_BUILD_TOOLS=ON
cmake --build cmake-build-debug --target RBFs sod_validation -j"$(nproc)"
./cmake-build-debug/RBFs
./cmake-build-debug/sod_validation --cells 50,100,200,400 \
    --out results/sod_validation
```
