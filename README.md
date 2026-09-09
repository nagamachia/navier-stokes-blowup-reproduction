# Navier–Stokes Blow-up Reproduction

A low-cost numerical exploration of the finite-time singularity construction announced by OpenAI for the 3D incompressible Navier–Stokes equations.

> [!IMPORTANT]
> This repository is an independent numerical exploration. It does **not** claim to prove, formally verify, or independently validate the announced mathematical result. The goal is to reproduce and study numerically accessible pre-singular scaling and cancellation mechanisms.

## Goals

1. Extract the numerically relevant construction from the published proof.
2. Reproduce the self-similar core scaling before building a full CFD solver.
3. Build and validate a small 3D incompressible pseudo-spectral solver.
4. Compare 64^3, 128^3, and 256^3 runs on a personal computer.
5. Measure maximum velocity, maximum vorticity, energy, enstrophy, dissipation, and residuals.
6. Investigate whether dynamic rescaling can follow the shrinking core more efficiently than brute-force DNS.

## Constraints

The project deliberately targets inexpensive personal computing:

- CPU first; no GPU requirement
- C++20
- FP64
- OpenMP for shared-memory parallelism
- FFTW for practical spectral transforms
- MPI only if a later experiment genuinely requires multiple nodes

## Roadmap

### Phase 0 — Paper-to-code

Extract definitions, coordinates, scaling laws, assumptions, and observables from the published construction. Implement the smallest possible model of the self-similar core without solving the full Navier–Stokes PDE.

### Phase 1 — Numerical infrastructure

Add tests, CSV output, reproducible parameter files, and basic convergence analysis.

### Phase 2 — Incompressible solver

Implement a periodic 3D Fourier pseudo-spectral Navier–Stokes solver with 2/3 dealiasing and validate it against standard test cases such as Taylor–Green vortex.

### Phase 3 — Blow-up construction

Introduce progressively more of the published construction. Compare theoretical scaling with numerical measurements and explicitly measure the cancellation/residual terms.

### Phase 4 — Dynamic rescaling

If direct resolution becomes the limiting factor, transform to coordinates that keep the shrinking core at approximately fixed computational size.

## Phase 0 quick start

The core-scaling executable checks only the asymptotic scaling algebra. It is intentionally **not** a Navier–Stokes solver.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/core_scaling
```

Optional arguments are

```bash
./build/core_scaling <h> <decades> <samples_per_decade> <output.csv>
```

For example:

```bash
./build/core_scaling 0.005 8 20 core_scaling.csv
```

The default experiment checks that the shrinking-core estimate gives

`U ~ tau^(-0.505)`

while

`E_core ~ tau^(0.485)`.

The velocity scale therefore grows as `tau -> 0`, while this core-energy estimate decreases.

## Solver verification

CTest covers the verification sequence in `docs/roadmap.md`:

1. Fourier differentiation on a periodic trigonometric function.
2. Spectral divergence-free projection, including idempotence and non-increasing modal energy.
3. 2/3 dealiasing using a deliberately aliased quadratic product.
4. Viscous decay of an incompressible Fourier mode against its analytic decay law.
5. A small 3D direct-DFT pseudo-spectral Taylor–Green regression using RK4, projection, rotational-form nonlinearity, and 2/3 dealiasing.

The direct DFT implementation is intentionally tiny and auditable. A separate FFTW implementation carries the same conventions into practical grid sizes and is also covered by a CI smoke test.

Run all checks with

```bash
ctest --test-dir build --output-on-failure
```

## FFTW resolution-study runner

Install FFTW development headers before configuring on Debian/Ubuntu:

```bash
sudo apt-get install libfftw3-dev
```

Then a first 64^3 Taylor–Green baseline can be run with

```bash
./build/fftw_taylor_green 64 0.1 0.001 0.1 taylor_green_64.csv
```

The CSV records maximum velocity, maximum vorticity, kinetic energy, enstrophy, viscous dissipation, divergence error, a discrete PDE consistency residual, projected nonlinear magnitude, timestep, and CFL number. The full Step 4 protocol and the criteria for moving from 64^3 to 128^3 and 256^3 are documented in `docs/resolution-study.md`.

## What counts as success?

The first milestone is deliberately modest: reproduce predicted pre-singular scaling over an increasing time interval as numerical resolution increases. A finite-resolution computation cannot demonstrate an actual infinite velocity at the singular time.

## Upstream material

- OpenAI announcement: https://openai.com/index/navier-stokes-solution/
- OpenAI paper: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- OpenAI formalization repository: https://github.com/openai/NavierStokesAndEuler

Theory/source notes are maintained in `docs/theory.md`. Exact similarity profiles and forcing will not be implemented until their definitions are pinned to the paper rather than inferred from secondary summaries.

## Status

- **Step 1:** initial auditable theory extraction complete; exact PDF equation-number pinning remains before full-profile implementation.
- **Step 2:** minimal core-scaling executable and release-active regression check complete.
- **Step 3:** all five solver-verification checks in `docs/roadmap.md` are implemented and passing in CI using an auditable direct-DFT reference implementation.
- **Step 4:** FFTW-backed 3D runner, required diagnostics, CI smoke regression, and the resolution-study protocol are in place. The next numerical milestone is a controlled 64^3 baseline and timestep-convergence study before considering 128^3.
