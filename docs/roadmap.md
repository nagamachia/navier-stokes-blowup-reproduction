# Experiment Roadmap

## Design principle

Spend mathematical effort before compute. We should not increase grid size until the previous level has a measurable failure caused by resolution rather than by an implementation error.

## Step 1 — Auditable theory extraction

Deliverable: `docs/theory.md`

For every quantity later implemented, record:

- notation used in the paper;
- section/equation reference;
- dimensional or nondimensional interpretation;
- asymptotic scaling as the singular time is approached;
- whether the expression is exact, asymptotic, or our numerical simplification.

No guessed profile should be labeled as the OpenAI construction.

## Step 2 — Minimal core experiment

Deliverable: `src/core_scaling.cpp`

Start without FFTW and without solving Navier–Stokes. Evaluate the extracted similarity/scaling model on a small grid and output CSV diagnostics. This separates understanding of the construction from PDE-solver errors.

Initial target: seconds to minutes on an ordinary CPU.

## Step 3 — Solver verification

Deliverable: periodic pseudo-spectral solver plus tests.

Order:

1. Fourier differentiation tests.
2. Divergence-free projection test.
3. 2/3 dealiasing test.
4. Viscous decay of a known Fourier mode.
5. Taylor–Green vortex regression test.

Only after these pass should the blow-up construction be inserted.

## Step 4 — Resolution study

Nominal sequence:

- 64^3 for debugging;
- 128^3 for routine experiments;
- 256^3 only when justified by convergence data.

Track at minimum:

- max velocity;
- max vorticity;
- kinetic energy;
- enstrophy;
- viscous dissipation;
- divergence error;
- PDE residual and relevant term-by-term cancellations;
- timestep and CFL number.

## Step 5 — Dynamic rescaling

If the physical core becomes under-resolved, test a dynamically rescaled formulation before paying for substantially larger DNS. The aim is to keep the shrinking structure O(1) in computational coordinates.

## Compute budget policy

Default budget: existing personal computer and open-source software. GPU, cloud HPC, and MPI are deferred until measurements show that they solve a specific bottleneck.

## Scientific caution

A numerical run can provide evidence about pre-singular behavior and consistency with a scaling law. It cannot by itself establish a mathematical finite-time singularity. Conversely, failure to reproduce a singular regime at low resolution is not evidence that the analytical construction is false.
