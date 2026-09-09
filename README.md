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
- FFTW when the spectral solver is introduced
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

## What counts as success?

The first milestone is deliberately modest: reproduce predicted pre-singular scaling over an increasing time interval as numerical resolution increases. A finite-resolution computation cannot demonstrate an actual infinite velocity at the singular time.

## Upstream material

- OpenAI announcement: https://openai.com/index/navier-stokes-solution/
- OpenAI formalization repository: https://github.com/openai/NavierStokesAndEuler

The exact paper equations used by code in this repository will be recorded with section/equation references in `docs/theory.md` before implementation.

## Status

**Step 1 in progress:** extract a minimal, auditable paper-to-code specification.
