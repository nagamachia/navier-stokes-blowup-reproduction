# Resolution study protocol

This document implements Step 4 of `docs/roadmap.md` as a reproducible numerical protocol. It deliberately starts with the Taylor–Green regression problem before any blow-up profile is inserted.

## Purpose

The resolution study answers a narrow question: is a discrepancy caused by spatial resolution, or is it still explainable by implementation or timestep error?

Grid size must not be increased merely because a larger run is available. The default sequence is:

1. `64^3` for debugging and the first convergence baseline;
2. `128^3` only after the 64^3 run is stable and its dominant error is plausibly spatial;
3. `256^3` only after 64^3 versus 128^3 data shows a measurable resolution-dependent discrepancy worth resolving.

The direct-DFT tests remain the auditable reference for spectral conventions. The FFTW runner is the practical implementation for larger grids.

## Build

On Debian/Ubuntu, install FFTW development headers first:

```bash
sudo apt-get install libfftw3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

All regression tests should pass before a resolution run is interpreted.

## FFTW Taylor–Green runner

The executable is:

```bash
./build/fftw_taylor_green <N> <final_time> <dt> <viscosity> <output.csv>
```

A first 64^3 baseline is:

```bash
./build/fftw_taylor_green 64 0.1 0.001 0.1 taylor_green_64.csv
```

If that run is stable and timestep error has been checked, the corresponding 128^3 run is:

```bash
./build/fftw_taylor_green 128 0.1 0.0005 0.1 taylor_green_128.csv
```

Do not infer spatial convergence from these two commands alone, because the timestep also changes. For a clean spatial comparison, first run both resolutions with a timestep small enough that halving `dt` at fixed `N` changes the observables much less than changing `N`.

## Numerical formulation currently exercised

The FFTW runner uses:

- a periodic `2*pi` box in all three directions;
- complex-to-complex FFTW transforms with the forward transform normalized by the number of grid points;
- spectral derivatives using integer Fourier wave numbers;
- incompressibility projection `P_k = I - k k^T / |k|^2` for nonzero modes;
- component-wise 2/3 dealiasing cutoff;
- rotational-form nonlinearity `u x omega`, followed by projection;
- explicit fourth-order Runge–Kutta time integration;
- spectral viscosity `-nu |k|^2 u_hat`.

The current Taylor–Green initial condition is a two-dimensional Taylor–Green vortex embedded in the 3D box. For this field, the projected nonlinear term vanishes and the exact velocity amplitude decays as `exp(-2 nu t)`. This makes it useful as an end-to-end regression target, but it is not yet a demanding turbulent resolution benchmark.

## CSV diagnostics

`fftw_taylor_green` writes one row per timestep with these columns:

- `step`: integer timestep index;
- `time`: physical simulation time;
- `max_velocity`: maximum pointwise velocity magnitude;
- `max_vorticity`: maximum pointwise vorticity magnitude;
- `kinetic_energy`: spectral kinetic energy `0.5 * sum |u_hat|^2`;
- `enstrophy`: spectral enstrophy `0.5 * sum |omega_hat|^2`;
- `viscous_dissipation`: `2 * nu * enstrophy`;
- `divergence_l2`: spectral L2 norm of `k dot u_hat`;
- `pde_residual_l2`: discrete consistency residual `(u_n-u_{n-1})/dt - RHS(u_n)`; this is a timestep diagnostic, not an exact continuous PDE residual;
- `projected_nonlinear_l2`: L2 norm of the projected, dealiased rotational nonlinear term before viscosity;
- `dt`: timestep;
- `cfl`: `max_velocity * dt / dx`, with `dx = 2*pi/N`.

These cover the minimum observables required by Step 4, with one caveat: term-by-term cancellation diagnostics for the announced blow-up construction cannot be meaningful until the exact construction and forcing are implemented from the paper.

## Acceptance checks before increasing N

For each resolution, perform at least one timestep-halving run. Increasing `N` is justified only after all of the following are true:

1. all CTest regression tests pass;
2. `divergence_l2` remains near floating-point noise relative to the velocity spectrum;
3. CFL remains comfortably below the explicit-stability limit used for the experiment;
4. halving `dt` changes the observables substantially less than increasing `N`;
5. no diagnostic grows solely because of an obvious normalization or aliasing error;
6. the quantity motivating the larger grid shows a systematic resolution dependence.

For the current Taylor–Green exact solution, also compare the final numerical velocity against `exp(-2 nu t)` times the initial field. The executable exits nonzero if its built-in exact-solution or divergence regression exceeds tolerance.

## What Step 4 does not yet establish

A converged Taylor–Green run validates the numerical infrastructure, not the finite-time singularity construction. Before the blow-up construction is inserted, `docs/theory.md` must be extended with exact paper section/equation references for every implemented profile, coordinate transform, pressure term, oscillatory correction, and forcing term.

Once those expressions are pinned to the primary source, the same diagnostics and resolution protocol can be reused for the pre-singular experiment.
