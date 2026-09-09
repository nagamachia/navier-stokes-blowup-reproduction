# Theory notes for the numerical reproduction

This document separates three levels of confidence:

- **Official result statement**: checked against OpenAI's announcement and formalization repository.
- **Published construction scaling**: used only where the scaling law is consistently reported from the paper; exact equation numbers are still to be pinned to the PDF before any claim of exact reproduction.
- **Numerical toy model**: deliberately simpler than the published construction and never labeled as the OpenAI solution itself.

## 1. Problem being explored

The announced Navier–Stokes result concerns the 3D incompressible equations with positive viscosity and smooth forcing. OpenAI states that there are smooth data and forcing for which a finite-time singularity develops. The official Lean repository states two versions: one on the whole space R^3 and one on the periodic torus R^3/Z^3.

For our first numerical work we do **not** attempt to reproduce the full proof construction. The immediate target is the pre-singular scaling of the concentrating core.

## 2. Normalized singular time

Use the normalized singular time

`T = 1`

and define remaining time

`tau = 1 - t`.

The numerical experiments always stop at `tau > 0`. No finite-resolution run will evaluate the singular time itself.

## 3. Core scaling used in the first experiment

The concentrating-core scaling reported from the published construction is

`L_r(tau) ~ tau^(1/2)`

for the two transverse/radial directions and

`L_z(tau) ~ tau^(1/2 - h)`

for the axial direction, with

`0 < h < 1/100`.

The characteristic large velocity components scale as

`U(tau) ~ tau^(-1/2 - h)`.

The radial velocity is reported at the weaker scale

`U_r(tau) = O(tau^(-1/2))`.

These are asymptotic scaling laws, not yet a complete velocity profile.

## 4. Energy consistency check

A core with two radial dimensions of scale `L_r`, one axial dimension of scale `L_z`, and velocity amplitude `U` has the order-of-magnitude energy

`E_core ~ U^2 * L_r^2 * L_z`.

Substituting the scalings above gives

`E_core ~ tau^(1/2 - 3h)`.

Since `h < 1/100`, the exponent is positive, so this core-energy estimate tends to zero even while the pointwise velocity amplitude diverges. This is the first scaling relation we can test without a Navier–Stokes solver.

## 5. Why this is not yet a reproduction of the proof

A shrinking singular background field by itself is not enough. When inserted into Navier–Stokes, the momentum residual would generally become singular. OpenAI's announcement emphasizes that the acceleration, pressure-gradient, nonlinear transport, and viscous terms become large but cancel in a precise way so that the applied force remains smooth.

Public descriptions of the paper report that the construction adds oscillatory corrections/pulses to cancel the singular residual. We will not implement these corrections until their exact definitions have been extracted directly from the paper.

Therefore `src/core_scaling.cpp` will only test scaling algebra. It does not solve the PDE and does not model the residual-cancellation mechanism.

## 6. Quantities for Phase 0

For logarithmically decreasing `tau`, record

- `tau`
- `t = 1 - tau`
- `L_r = tau^(1/2)`
- `L_z = tau^(1/2-h)`
- `U = tau^(-1/2-h)`
- `U_r_scale = tau^(-1/2)`
- `V_core = L_r^2 * L_z`
- `E_core = U^2 * V_core`
- analytic reference `tau^(1/2-3h)`
- relative error between the algebraic energy expression and its analytic power law

This first executable should reproduce the exponent exactly up to floating-point roundoff.

## 7. Parameter choice for the toy experiment

Use

`h = 0.005`

by default. This lies safely inside `0 < h < 0.01` and yields

`U ~ tau^(-0.505)`

and

`E_core ~ tau^(0.485)`.

The choice is for visualization and numerical sanity checks only; it is not claimed to be the parameter used in the exact construction.

## 8. Source status

Primary sources:

- OpenAI announcement: https://openai.com/index/navier-stokes-solution/
- OpenAI formalization repository: https://github.com/openai/NavierStokesAndEuler
- Paper URL published by OpenAI: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf

Before implementing the actual similarity profiles or forcing, this document must be extended with exact PDF section/equation references for every expression copied into code.

## 9. Next paper-extraction targets

The next exact-extraction pass should identify, directly from the paper:

1. the similarity coordinates;
2. the definitions of the background/core velocity components;
3. the profile functions used in the core;
4. the pressure definition;
5. the annular transition region;
6. the oscillatory pulse families;
7. the residual/stress cancellation identities;
8. the precise smooth forcing constructed at the end.

Until those are pinned down, Phase 0 remains a scaling experiment rather than a full reproduction.
