# 数値再現のための理論ノート

この文書では数式の確度を次の3段階に分ける。

- **論文に明記された厳密な式**: OpenAI 論文から節・式番号付きで確認したもの
- **数値実装上の帰結**: 上記から代数的・漸近的に導くもの
- **数値トイモデル**: 論文構成より意図的に単純化したもの。OpenAI の exact construction とは呼ばない

原典: `FINITE TIME BLOWUP FOR NAVIER–STOKES`, OpenAI, 2026。

## 1. 問題設定と規格化

Theorem 1.1 と式 (1.1) では、任意の正粘性に対して、滑らかでコンパクト台を持つ forcing と、静止状態から始まり、`t < 1` では運動エネルギーが一様有界である一方、`t -> 1-` で速度が非有界になる 3D 非圧縮 Navier–Stokes 解を構成する。

粘性1の構成では

`tau = 1 - t`

と置く。本リポジトリの数値計算は常に `tau > 0` で停止する。

## 2. 厳密な類似座標

Section 3.1 と Section 4.1 で

`A = 1/2 + h`

`D = 1/2 - h`

を用いる。最終構成では `0 < h < 1/100` で、`A + D = 1`。

式 (3.2)、(4.1):

`z = q^D eta`

`tau = q (1 - eta^2)`

`d = 1 - eta^2`

`L = 1 - 2 h eta^2`

`s = r^2 / 2`

`X = s/q = r^2/(2q)`

ここで `q > 0`, `|eta| < 1`。

`eta = z q^(-D)` を消去すると

`q - z^2 q^(2h) = tau`

を得る。これは Section 4.1 の本文にも明記されている。物理 branch では

`d/dq [q - z^2 q^(2h)] = 1 - 2 h z^2 q^(2h-1) = L >= 1 - 2h > 0`

となるため、正解は一意。

> [!NOTE]
> 初期の転記で `q - z^2/q^(2h) = tau` としていたが、Lemma 4.1 の有限差分照合で不整合を検出し、PDF の Section 4.1 を再確認して修正した。

実装:

- `include/similarity_coordinates.hpp`
- `src/similarity_coordinates.cpp`
- `tests/similarity_coordinates.cpp`

## 3. Lemma 4.1 の微分作用素

式 (4.2) は滑らかな profile `f(X,eta)` と実数 `b` に対して

`partial_t(q^b f) = q^(b-1) T_b f`

`partial_z(q^b f) = q^(b-D) Z_b f`

を与える。

`D_X f = X partial_X f`

として

`T_b f = L^(-1) (-b f + D eta partial_eta f + D_X f)`

`Z_b f = L^(-1) (2 b eta f + d partial_eta f - 2 eta D_X f)`

である。

論文の proof では

`q_t = -L^(-1)`

`eta_t = D eta/(qL)`

`X_t = X/(qL)`

`q_z = 2 eta q^(1-D)/L`

`eta_z = d/(q^D L)`

`X_z = -2 eta X/(q^D L)`

を導いている。

実装:

- `include/similarity_operators.hpp`
- `tests/similarity_operators.cpp`

テストでは manufactured profile を使い、物理座標で直接評価した中心有限差分 `partial_t`, `partial_z` と `T_b`, `Z_b` を比較する。

## 4. leading velocity と pressure

式 (4.3):

`u_theta^(0) = q^(-A) E(X,eta)`

`u_z^(0) = q^(-A) U(X,eta)`

`r u_r^(0) = V0(X,eta)`

`p^(0) = q^(-2A) Pi(X,eta)`

式 (4.4) の軸上正則性:

`E = sqrt(2X) F`

`V0 = X v0`

で、`F, U, v0, Pi` は inner profile rectangle 上で滑らか。論文ではさらに `E = C^(-1) sqrt(2X) phi` と因子分解する。

式 (4.5) は Cartesian 成分を与える。`r=0` 近傍の数値実装では円柱基底を直接扱うより、この Cartesian 表現を使う方が安全。

## 5. 非圧縮条件と圧力

式 (4.6) の radial average:

`A_X(f)(X,eta) = (1/X) integral_0^X f(x,eta) dx`

軸上では

`A_X(f)(0,eta) = f(0,eta)`

と滑らかに延長する。

式 (4.7):

`V0 = (X/L) [2 eta U - 2 D eta A_X(U) - d partial_eta A_X(U)]`

および

`partial_X Pi = E^2/(2X)`

を与える。

Section 3.1 の pressure normalization は

`Pi(X,eta) = - integral_X^infinity E(x,eta)^2/(2x) dx`

である。

次の直接実装対象は、任意の tabulated `E,U` に対して `A_X`, `V0`, `Pi` を再構成し、manufactured profile で積分精度と differential identity を検証すること。

## 6. コア形状と blow-up rate

固定された inner profile 領域で `q ~ tau` なので

`ell_r ~ tau^(1/2)`

`ell_z ~ tau^(1/2-h)`

となる。

主要成分は

`||u_theta^(0)||_inf ~ tau^(-1/2-h)`

`||u_z^(0)||_inf ~ tau^(-1/2-h)`

一方

`||u_r^(0)||_inf = O(tau^(-1/2))`

である。

固定 `0 < X_* < X_c` に対して

`z = 0`

`r = sqrt(2 X_* tau)`

とすると `q=tau`, `eta=0` なので

`u_theta^(0) = E(X_*,0) tau^(-1/2-h)`

となる。

Theorem 3.1(iv)、Proposition 9.9 では全補正後にも

`u_theta,loc(sqrt(2 X_in tau),0,0,1-tau) = tau^(-A) [e0 + O(tau^(2h))]`

という leading growth が保たれる。exact construction 実装後の主要回帰観測量にする。

## 7. エネルギー整合性

コア体積は

`tau^(3/2-h)`

主要速度スケールは

`tau^(-1/2-h)`

なので、コア運動エネルギーは

`E_core ~ tau^(1/2-3h)`

となる。`h < 1/100` では指数が正なので、点wise速度が増大しながらコアのエネルギー寄与は0へ向かい得る。

`src/core_scaling.cpp` はこの代数スケーリングのみを検証するトイモデル。

## 8. leading profile は単純な閉形式ではない

論文は `E(X,eta), U(X,eta)` を単純な Gaussian や短い解析式として与えていない。

Section 4 の Theorem 4.6 は、適切な

`h in (0,1/100)`

`lambda > 0`

`C > 1`

`0 < X_a < X_b`

と、regular-axis、annular-stress、matching、moment、heat-exterior 条件を満たす profile を構成する。

Appendix B が inner profile、Appendix A が outer profile と heat exterior の構成・matching を担う。

したがって exact numerical reproduction では、便利な surrogate profile を OpenAI 構成として扱わず、この constructive procedure 自体を段階的に数値化する。

## 9. annular residual と oscillatory cancellation

inner field と exterior の接続により active annulus に singular momentum residual が残る。

active annulus は類似座標で概ね

`X_a < X < X_b`

に固定される。Theorem 4.6 と式 (4.24) を参照。

式 (4.7)–(4.11) は leading stress quantities を定義し、式 (4.23) は admissible stress cone 条件を与える。

oscillatory pulses の averaged quadratic momentum flux が不足 stress を実現する。

- Proposition 7.5: leading stress の wave covariance realization
- Proposition 7.6: higher-order stress correction
- Proposition 9.6: residual improvement
- Proposition 9.9: summed local field

実装順序は

1. leading profile と stress
2. wave realization
3. residual-improvement correction

とする。

## 10. exterior と最終 smooth forcing

active radial profile の外側では純方位角 heat-flow exterior を保つ。Section 3.1 と式 (4.29) を参照。

Section 10 では local construction を物理空間・時間で局所化する。

- Proposition 10.1: vector potential レベルで spatial cutoff を行い非圧縮性を保存
- temporal cutoff により初期速度を0にする
- 最終 force は局所化された `(u,p)` の Navier–Stokes momentum residual として定義
- Lemma 10.3: `t=1` をまたぐ smooth extension と compact support

faithful reproduction では forcing を独立に仮定せず、完成した `(u,p)` から residual として計算する。

一般粘性では

`u_nu(x,t) = sqrt(nu) u(x/sqrt(nu),t)`

`p_nu(x,t) = nu p(x/sqrt(nu),t)`

`f_nu(x,t) = sqrt(nu) f(x/sqrt(nu),t)`

を使う。

## 11. 数値回帰の階層

1. `q(z,tau)` を解き式 (3.2)/(4.1) を再構成 — **実装済み**
2. `T_b`, `Z_b` を物理座標有限差分と比較 — **実装済み**
3. tabulated `E,U` から `A_X`, `V0`, `Pi` を計算
4. 式 (4.4)–(4.5) の軸上正則性を確認
5. `r=sqrt(2 X_in tau), z=0` 経路で指数 `-A` を fit
6. inner / active annulus / heat exterior ごとに momentum residual を評価
7. waves 実装後、background residual と averaged quadratic wave stress のキャンセルを測定
8. `tau -> 0` で最終 residual とその微分の smooth/flat behavior を確認

## 12. 現在の実装境界

推測なしで実装可能なところまで式番号を固定済み:

- 類似座標: (3.2), (4.1) — **実装済み**
- 微分作用素: Lemma 4.1, (4.2) — **実装済み**
- leading velocity / pressure ansatz: (4.3)–(4.5)
- radial averaging、非圧縮条件、pressure balance: (4.6)–(4.7) — **次の実装対象**
- active annulus と leading-profile 性質: (4.24), Theorem 4.6
- blow-up sampling path: Theorem 3.1(iv), Proposition 9.9
- localization / final forcing: Proposition 10.1, Lemma 10.3

次は `A_X`, `V0`, `Pi` の汎用数値部品と manufactured-profile 回帰を実装する。その後、式 (4.8)–(4.11) と Sections 6–7 の pulse / stress realization を追加抽出する。
