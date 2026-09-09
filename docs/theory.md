# Theory notes for the numerical reproduction

This document separates three levels of confidence:

- **Paper-exact identities**: transcribed from the OpenAI paper with section/equation references.
- **Consequences for numerical work**: algebraic or asymptotic consequences of those identities.
- **Numerical toy models**: deliberately simpler than the published construction and never labeled as the OpenAI solution itself.

Primary paper: `FINITE TIME BLOWUP FOR NAVIER–STOKES`, OpenAI, 2026.

## 1. Problem and normalization

The paper constructs, for every positive viscosity, a smooth compactly supported force and a solution of the 3D incompressible Navier–Stokes equations that starts from rest, has uniformly bounded kinetic energy for `t < 1`, and develops unbounded velocity as `t -> 1-`; see Theorem 1.1 and equation (1.1).

For the construction at viscosity one, define

`tau = 1 - t`.

The proof outline rescales the viscosity-one construction to arbitrary `nu > 0`; see Section 3 and equations (10.22)-(10.23). Numerical experiments in this repository stop at `tau > 0`.

## 2. Exact similarity coordinates

The concentrating leading field is introduced in Section 3.1 and developed in Section 4.1.

The fixed exponents are

`A = 1/2 + h`,

`D = 1/2 - h`,

with the final construction choosing `0 < h < 1/100`. Hence `A + D = 1`.

Equation (3.2), repeated as part of (4.1), defines the similarity coordinates by

`tau = q (1 - eta^2)`,

`z = q^D eta`,

`X = r^2 / (2 q)`,

with `q > 0` and `-1 < eta < 1`.

Equivalently, eliminating `eta`, `q` is the unique solution of

`q - z^2 / q^(2h) = tau`.

On fixed compact subsets of similarity coordinates away from `eta = +/-1`, `q ~ tau`. For bounded `X`, the limit `q -> 0` approaches the singular point.

Section 4.1 also introduces

`d = 1 - eta^2`,

`L = 1 - 2 h eta^2`,

`s = r^2 / 2`,

so `X = s/q`; see equation (4.1).

## 3. Exact derivative operators

Lemma 4.1, equation (4.2), gives the coordinate calculus. For a smooth profile `f(X,eta)` and real `b`,

`partial_t(q^b f) = q^(b-1) T_b f`,

`partial_z(q^b f) = q^(b-D) Z_b f`,

where

`T_b f = L^(-1) (-b f + D eta partial_eta f + D_X f)`,

`Z_b f = L^(-1) (2 b eta f + d partial_eta f - 2 eta D_X f)`,

and

`D_X f = X partial_X f`.

These formulas should be implemented and unit-tested before a numerical implementation of the exact leading profile.

## 4. Leading velocity and pressure fields

The exact leading-field ansatz is equation (4.3):

`u_theta^(0) = q^(-A) E(X,eta)`,

`u_z^(0) = q^(-A) U(X,eta)`,

`r u_r^(0) = V0(X,eta)`,

`p^(0) = q^(-2A) Pi(X,eta)`.

For regularity at the cylindrical axis the azimuthal profile is factored as

`E = C^(-1) sqrt(2X) phi`,

with fixed `C > 1`. Equation (4.4) states the axis regularity structure

`E = sqrt(2X) F`,

`V0 = X v0`,

with `F = phi/C` and `F, U, v0, Pi` smooth on the closed inner profile rectangle.

Equation (4.5) gives the corresponding Cartesian components. This is the preferred formula for code evaluating the field near `r = 0`, because it avoids treating the cylindrical basis as defined on the axis.

## 5. Incompressibility and pressure identities

The profiles `E` and `U` are the primary choices. Incompressibility determines `V0`, while radial centrifugal balance determines the radial pressure derivative.

Define the radial average, equation (4.6),

`A_X(f)(X,eta) = (1/X) integral_0^X f(x,eta) dx`,

with its smooth axis value `A_X(f)(0,eta) = f(0,eta)`.

Then equation (4.7) gives

`V0 = (X/L) [2 eta U - 2 D eta A_X(U) - d partial_eta A_X(U)]`,

and

`partial_X Pi = E^2 / (2X)`.

The pressure normalization used in the proof is equivalently

`Pi(X,eta) = - integral_X^infinity E(x,eta)^2 / (2x) dx`;

see Section 3.1 immediately after (3.2).

These identities are direct implementation targets: a discretized profile should satisfy incompressibility and radial pressure balance to numerical tolerance before the full residual is evaluated.

## 6. Core geometry and blow-up rate

For fixed inner profile bounds `0 <= X <= X_c` and `|eta| <= eta_c < 1`, Section 3.1 defines the physical core `C_tau`. Since `q ~ tau` there,

`ell_r ~ tau^(1/2)`,

`ell_z ~ tau^(1/2-h)`.

The leading tangential components obey

`||u_theta^(0)||_inf ~ tau^(-1/2-h)`,

`||u_z^(0)||_inf ~ tau^(-1/2-h)`,

while

`||u_r^(0)||_inf = O(tau^(-1/2))`.

A particularly useful exact numerical sampling path is already identified by the proof. For any fixed `0 < X_* < X_c`, take

`z = 0`, `r = sqrt(2 X_* tau)`.

Then `q = tau`, `eta = 0`, and

`u_theta^(0) = E(X_*,0) tau^(-1/2-h)`.

After all corrections are summed, Theorem 3.1(iv), proved in Proposition 9.9, preserves this leading growth in the form

`u_theta,loc(sqrt(2 X_in tau), 0, 0, 1-tau) = tau^(-A) [e0 + O(tau^(2h))]`.

This path should be the primary regression observable once the construction beyond the toy scaling model is implemented.

## 7. Energy consistency

The core volume has order

`tau^(3/2-h)`.

Combining this with the dominant velocity scale `tau^(-1/2-h)` gives the core kinetic-energy order

`E_core ~ tau^(1/2-3h)`.

Because the construction has `h < 1/100`, this exponent is positive. Thus the core contribution can tend to zero while the pointwise velocity diverges. This agrees with the physical description in Section 2.1 and Section 3.1.

The existing `src/core_scaling.cpp` tests only this scaling algebra; it is not the leading-profile PDE construction.

## 8. Leading profile is not an arbitrary analytic formula

The paper does not simply provide a short closed-form pair `E(X,eta), U(X,eta)` that can be copied into a CFD initial condition. The profiles are constructed through matching, moment constraints, an inner analytic construction, a heat exterior, and stress-cone inequalities.

Section 4.4 culminates in Theorem 4.6. It asserts existence of fixed

`h in (0,1/100)`, `lambda > 0`, `C > 1`, `0 < X_a < X_b`,

and profiles `E, U, Pi` satisfying the required regular-axis, annular-stress, matching, moment, and exterior properties. Appendix B constructs the inner profile; Appendix A constructs/matches the outer profile and heat exterior.

Therefore an exact numerical reproduction must implement the constructive profile procedure, not invent a convenient Gaussian or Taylor-Green surrogate and call it the paper's profile.

## 9. Annular residual and oscillatory cancellation

The background profile alone is insufficient. The physical description in Section 2 and proof outline in Section 3 explain that joining the concentrating inner field to the exterior leaves a singular momentum residual in an annulus.

The active annulus is fixed in similarity coordinates:

`X_a < X < X_b`, `-1 <= eta <= 1`;

see Theorem 4.6 and equation (4.24) for the inner/outer collar definitions.

The construction rewrites the tangential residual as a radial stress divergence. Equations (4.7)-(4.11) define the leading stress quantities, and the admissible stress-cone condition is expressed in equation (4.23). Theorem 4.6 arranges the leading profile so this stress is realizable by the oscillatory wave families.

The oscillatory pulses are not optional numerical decoration: their averaged quadratic momentum flux supplies the missing stress. Proposition 7.5 realizes the leading stress by positive wave covariances; Proposition 7.6 supplies higher-order stress corrections. Subsequent corrections improve the remaining residual, culminating in Proposition 9.6 and the summed local field of Proposition 9.9.

This establishes the implementation order: leading profile and stress first, then wave realization, then residual-improvement corrections.

## 10. Exterior and final smooth forcing

Outside the active radial profile region, the construction preserves an exact purely azimuthal heat-flow exterior; Section 3.1 points to equation (4.29). Its momentum residual vanishes there.

After the local construction is summed, Section 10 localizes it in physical space and time. Proposition 10.1 applies the spatial cutoff at the vector-potential level so incompressibility is preserved, and applies a temporal cutoff so the initial velocity is zero.

The final force is the Navier-Stokes momentum residual of the localized fields. Lemma 10.3 supplies the smooth extension through `t = 1` with compact space-time support. Thus the force should not be approximated independently: in a faithful reproduction it is computed from the completed localized `(u,p)` and checked for smooth extension.

For arbitrary viscosity, use the scaling stated in Section 3:

`u_nu(x,t) = sqrt(nu) u(x/sqrt(nu),t)`,

`p_nu(x,t) = nu p(x/sqrt(nu),t)`,

`f_nu(x,t) = sqrt(nu) f(x/sqrt(nu),t)`.

## 11. Numerically testable invariants and diagnostics

The exact extraction suggests the following regression hierarchy:

1. Solve the implicit similarity-coordinate equation for `q(z,tau)` and verify (3.2)/(4.1).
2. Unit-test the `T_b` and `Z_b` derivative identities against finite differences or automatic differentiation.
3. Given tabulated `E,U`, compute `V0` and `Pi` from (4.6)-(4.7) and verify Cartesian divergence and radial pressure balance.
4. Verify axis regularity using the factorizations (4.4)-(4.5).
5. Measure the leading path `r = sqrt(2 X_in tau), z=0` and fit the exponent `-A`.
6. Evaluate the momentum residual separately in the inner region, active annulus, and heat exterior.
7. Once waves are implemented, measure cancellation between background residual and averaged quadratic wave stress.
8. Track the final residual and its derivatives as `tau -> 0`; the target is smooth/flat behavior rather than merely small velocity error.

## 12. Implementation boundary after this extraction

The following items are now pinned directly to the paper and may be implemented without guessing:

- similarity coordinates: (3.2), (4.1);
- derivative operators: Lemma 4.1, (4.2);
- leading velocity/pressure ansatz: (4.3)-(4.5);
- radial averaging, incompressibility and pressure balance: (4.6)-(4.7);
- active annulus and leading-profile existence/properties: (4.24), Theorem 4.6;
- preserved blow-up sampling path: Theorem 3.1(iv), Proposition 9.9;
- localization/final forcing architecture: Proposition 10.1 and Lemma 10.3.

The next paper-extraction pass needed before coding the oscillatory realization should transcribe the exact definitions (4.8)-(4.11), the pulse coordinate/phase/amplitude formulas in Sections 6-7, and the stress realization of Propositions 7.5-7.6. Those formulas are substantially more involved and should be implemented only after their dependencies are pinned in the same way.
