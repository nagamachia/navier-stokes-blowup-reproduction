# Appendix A 数値化ノート

原典: OpenAI, *Finite Time Blowup for Navier–Stokes* (2026), Appendix A および Section 4.2。

この文書は Appendix A の constructive procedure を、推測した surrogate profile に置き換えず、実装可能な数値部品へ分解する。

## 1. five cumulative radial moments — 式 (4.15)

profile `E(X,eta), U(X,eta)` と angular momentum profile

`H = sqrt(2X) E`

に対して

`M(X,eta) = integral_0^X U dx`

`I(X,eta) = integral_0^X H dx`

`J(X,eta) = integral_0^X U H dx`

`S(X,eta) = integral_0^X (U^2 - E^2/2) dx`

`Cp(X,eta) = integral_0^X E^2/(2x) dx`

を使う。pressure は共通の axis datum `Pi_ax(eta)` に対して

`Pi = Pi_ax + Cp`

となる。

regular-axis factorization `E=sqrt(2X)F` を用いると

`H = 2 X F`

`E^2/(2X) = F^2`

となり、`Cp` の integrand は `X=0` でも正則になる。

実装: `include/leading_profile_relations.hpp::cumulative_radial_moments`。

## 2. profile join の意味 — Lemma 4.4

joining radius `Xh` より外で2つの profile `(U1,E1)`, `(U2,E2)` が一致し、`Xh` で上の5 moments が一致すれば、それより外側では pressure、`V0`、`Qs`、`Ns`、`ps`、shear、leading stress も一致する。

したがって Appendix A/B/C の局所編集では、編集後に5 momentsを元の値へ戻すことが exterior を壊さないための主要な数値拘束になる。

## 3. Lemma A.1 — moment matrix

互いに離れた正の radial intervals `Ij` に非負 smooth bump `beta_j` を置き、相異なる power `alpha_i` に対して

`B_ij = integral x^(alpha_i) beta_j(x) dx`

を作る。Lemma A.1 はこの行列の可逆性を与える。

実装:

- `include/moment_corrections.hpp::power_moment_matrix`
- partial-pivot Gaussian elimination による `solve_linear`
- 回帰テストでは3つの separated bump と3つの異なる power から非零 determinant と係数回収を確認する。

## 4. Lemma A.2 — quadratic moment correction

実際の moments は velocity profiles に二次依存するため、補正係数 `c` は

`F_eta(c) = B(eta)c + Q_eta(c,c) = d(eta)`

という有限次元系を満たす。原典が選ぶ zero-start branch に合わせて

`c_{n+1} = B^{-1}(d - Q(c_n,c_n)), c_0=0`

を実装する。

実装: `solve_small_quadratic_moment_system`。

現段階では各 `eta` 固定の有限次元問題を解く部品であり、`eta` に沿った `C^k` bounds の自動証明を行うものではない。今後は `eta` grid 上で係数関数を追跡し、有限差分で parameter derivative bounds を診断する。

## 5. Corollary A.3 への接続

原典では power-law correction patch 上で、2つの `U` bump と3つの `E` bump により5 momentsを調整する。Jacobian は axial 2×2 block と azimuthal 3×3 block に分かれ、異なる power weights により Lemma A.1 を適用する。

次の実装対象は以下。

1. `U0 = uc(eta)`, `E0 = e_* f_*(eta) x^alpha` の correction patch class
2. 5 bump coefficient から `(Delta M, Delta I, Delta J, Delta S, Delta Cp)` を返す exact quadratic map
3. numerical Jacobian と A.3 の block structure の照合
4. target discrepancy を与えた5係数 correction の回帰

これが通った後に Appendix A.2 の staged outer profile construction、heat exterior replacement、pressure datum の再現へ進む。

## 科学的注意

ここで実装しているのは論文に明記された moment matching mechanism である。一方、Appendix A.2 以降の profile schedule はまだ全て実装していないため、現時点のコードを Theorem 4.6 の完成 profile と呼ばない。
