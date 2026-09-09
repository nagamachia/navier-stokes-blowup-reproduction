# 数値再現のための理論ノート

この文書では、数式の確度を次の3段階に分ける。

- **論文に明記された厳密な式**: OpenAI 論文から節・式番号付きで転記したもの
- **数値実装上の帰結**: 上記から代数的・漸近的に導かれるもの
- **数値トイモデル**: 論文構成より意図的に単純化したもの。OpenAI の exact construction とは呼ばない

原典: `FINITE TIME BLOWUP FOR NAVIER–STOKES`, OpenAI, 2026。

## 1. 問題設定と規格化

Theorem 1.1 と式 (1.1) では、任意の正粘性に対して、滑らかでコンパクト台を持つ forcing と、静止状態から始まり、`t < 1` では運動エネルギーが一様有界である一方、`t -> 1-` で速度が非有界になる 3D 非圧縮 Navier–Stokes 解を構成する。

粘性1の構成では

`tau = 1 - t`

と置く。本リポジトリの数値計算では常に `tau > 0` で停止する。

一般の `nu > 0` へのスケーリングは Section 3 と式 (10.22)–(10.23) にある。

## 2. 厳密な類似座標

集中する leading field は Section 3.1 で導入され、Section 4.1 で展開される。

固定指数は

`A = 1/2 + h`

`D = 1/2 - h`

で、最終構成では `0 < h < 1/100`。したがって `A + D = 1`。

式 (3.2)、および式 (4.1) の一部として、類似座標は

`tau = q (1 - eta^2)`

`z = q^D eta`

`X = r^2 / (2q)`

で定義される。`q > 0`, `-1 < eta < 1`。

`eta` を消去すると `q` は

`q - z^2 / q^(2h) = tau`

の一意な正解である。

Section 4.1 ではさらに

`d = 1 - eta^2`

`L = 1 - 2 h eta^2`

`s = r^2 / 2`

`X = s/q`

を導入する。

この部分は `include/similarity_coordinates.hpp` に実装済みで、`tests/similarity_coordinates.cpp` に回帰テストがある。

## 3. 微分作用素

Lemma 4.1、式 (4.2) は、滑らかな `f(X,eta)` と実数 `b` に対して

`partial_t(q^b f) = q^(b-1) T_b f`

`partial_z(q^b f) = q^(b-D) Z_b f`

を与える。

ここで

`T_b f = L^(-1) (-b f + D eta partial_eta f + D_X f)`

`Z_b f = L^(-1) (2 b eta f + d partial_eta f - 2 eta D_X f)`

`D_X f = X partial_X f`

である。

exact leading profile を数値化する前に、この微分作用素を独立実装して有限差分などと照合する。

## 4. leading velocity と pressure

式 (4.3) の ansatz は

`u_theta^(0) = q^(-A) E(X,eta)`

`u_z^(0) = q^(-A) U(X,eta)`

`r u_r^(0) = V0(X,eta)`

`p^(0) = q^(-2A) Pi(X,eta)`

である。

円柱軸 `r=0` での正則性のため、方位角方向プロファイルは

`E = C^(-1) sqrt(2X) phi`

と因子分解される。式 (4.4) では

`E = sqrt(2X) F`

`V0 = X v0`

と書かれ、`F, U, v0, Pi` は inner profile rectangle 上で滑らか。

式 (4.5) は対応する Cartesian 成分を与える。`r=0` 近傍の実装では円柱基底を直接扱わず、式 (4.5) の Cartesian 表現を使う方が安全。

## 5. 非圧縮条件と圧力

`E` と `U` を主なプロファイルとして選ぶと、非圧縮条件が `V0` を決め、遠心力バランスが圧力の半径方向微分を決める。

式 (4.6) の半径平均:

`A_X(f)(X,eta) = (1/X) integral_0^X f(x,eta) dx`

軸上では滑らかな極限として

`A_X(f)(0,eta) = f(0,eta)`

を使う。

式 (4.7):

`V0 = (X/L) [2 eta U - 2 D eta A_X(U) - d partial_eta A_X(U)]`

および

`partial_X Pi = E^2 / (2X)`

を与える。

圧力の規格化は Section 3.1 の (3.2) 直後に

`Pi(X,eta) = - integral_X^infinity E(x,eta)^2 / (2x) dx`

として書かれている。

数値的には、tabulated `E,U` から `V0,Pi` を再構成し、Cartesian divergence と radial pressure balance を誤差評価するのが直接の実装対象。

## 6. コア形状と blow-up rate

固定された inner profile 領域

`0 <= X <= X_c`

`|eta| <= eta_c < 1`

では `q ~ tau`。

したがって物理長さは

`ell_r ~ tau^(1/2)`

`ell_z ~ tau^(1/2-h)`

となる。

leading tangential components は

`||u_theta^(0)||_inf ~ tau^(-1/2-h)`

`||u_z^(0)||_inf ~ tau^(-1/2-h)`

一方で

`||u_r^(0)||_inf = O(tau^(-1/2))`

である。

数値回帰に特に有用なのが、固定した `0 < X_* < X_c` に対する

`z = 0`

`r = sqrt(2 X_* tau)`

という経路。このとき

`q = tau`

`eta = 0`

なので

`u_theta^(0) = E(X_*,0) tau^(-1/2-h)`

となる。

補正を全て加えた後も Theorem 3.1(iv)、Proposition 9.9 により

`u_theta,loc(sqrt(2 X_in tau), 0, 0, 1-tau) = tau^(-A) [e0 + O(tau^(2h))]`

が保たれる。exact construction を実装した後の主要回帰観測量にする。

## 7. エネルギー整合性

コア体積は

`tau^(3/2-h)`

のオーダー。

主要速度スケール `tau^(-1/2-h)` を組み合わせると、コア運動エネルギーは

`E_core ~ tau^(1/2-3h)`

となる。

`h < 1/100` なので指数は正。したがって点wise速度が発散しても、コアのエネルギー寄与は0へ向かい得る。

`src/core_scaling.cpp` はこの代数スケーリングだけを検証するトイモデルであり、exact leading profile ではない。

## 8. leading profile は単純な閉形式ではない

論文は `E(X,eta), U(X,eta)` を単純な Gaussian や短い解析式として与えているわけではない。

Section 4.4 の Theorem 4.6 は、ある固定された

`h in (0,1/100)`

`lambda > 0`

`C > 1`

`0 < X_a < X_b`

と、必要な正則性、annular stress、matching、moment、exterior 条件を満たす `E,U,Pi` の存在を構成的に示す。

Appendix B が inner profile、Appendix A が outer profile と heat exterior の構成・matching を担う。

したがって exact numerical reproduction では、便利な surrogate profile を置くのではなく、この構成手順自体を数値化する必要がある。

## 9. annular residual と oscillatory cancellation

背景プロファイルだけでは不十分。Section 2 と Section 3 の説明では、集中する inner field と exterior を接続すると annulus に特異な momentum residual が残る。

active annulus は類似座標で

`X_a < X < X_b`

`-1 <= eta <= 1`

に固定される。Theorem 4.6 と式 (4.24) を参照。

式 (4.7)–(4.11) は leading stress quantities を定義し、式 (4.23) が admissible stress cone 条件を与える。

Theorem 4.6 は、この stress が後段の oscillatory wave families で実現可能になるよう leading profile を構成する。

oscillatory pulses は装飾ではなく、平均された二次運動量 flux が不足 stress を埋める本質的部分。

- Proposition 7.5: leading stress を正の wave covariance で実現
- Proposition 7.6: higher-order stress corrections
- Proposition 9.6: 残差改善
- Proposition 9.9: summed local field

したがって実装順序は

1. leading profile と stress
2. wave realization
3. residual-improvement corrections

となる。

## 10. exterior と最終 smooth forcing

active radial profile の外側では、exact な純方位角 heat-flow exterior が保たれる。Section 3.1 から式 (4.29) を参照。この領域では momentum residual は0。

Section 10 では local construction を物理空間・時間で局所化する。

- Proposition 10.1: vector potential レベルで空間 cut-off を行い、非圧縮性を保存
- 時間 cut-off により初期速度を0にする
- 最終 force は局所化された `(u,p)` の Navier–Stokes momentum residual として定義
- Lemma 10.3: `t=1` をまたぐ smooth extension と compact space-time support

したがって faithful reproduction では forcing を独立に近似するのではなく、完成した `(u,p)` から残差として計算する。

一般粘性では Section 3 の scaling

`u_nu(x,t) = sqrt(nu) u(x/sqrt(nu),t)`

`p_nu(x,t) = nu p(x/sqrt(nu),t)`

`f_nu(x,t) = sqrt(nu) f(x/sqrt(nu),t)`

を使う。

## 11. 数値的に検証可能な不変量・診断量

論文抽出から、次の回帰階層が得られる。

1. `q(z,tau)` を解き、式 (3.2)/(4.1) を再構成
2. `T_b`, `Z_b` を有限差分または自動微分と比較
3. tabulated `E,U` から `V0,Pi` を計算し、Cartesian divergence と radial pressure balance を検証
4. 式 (4.4)–(4.5) による軸上正則性を確認
5. `r = sqrt(2 X_in tau), z=0` 経路で指数 `-A` を fit
6. inner / active annulus / heat exterior ごとに momentum residual を評価
7. waves 実装後、background residual と averaged quadratic wave stress のキャンセルを測定
8. `tau -> 0` で最終 residual とその微分を追跡し、単に小さいだけでなく smooth/flat であることを確認

## 12. 現在の実装境界

論文から式番号まで固定できており、推測なしで実装してよい項目:

- 類似座標: (3.2), (4.1)
- 微分作用素: Lemma 4.1, (4.2)
- leading velocity / pressure ansatz: (4.3)–(4.5)
- radial averaging、非圧縮条件、pressure balance: (4.6)–(4.7)
- active annulus と leading-profile 性質: (4.24), Theorem 4.6
- blow-up sampling path: Theorem 3.1(iv), Proposition 9.9
- localization / final forcing: Proposition 10.1, Lemma 10.3

現在、類似座標は実装済み。

次に実装する順序:

1. Lemma 4.1 の `T_b`, `Z_b`
2. `A_X`、`V0`、`Pi`
3. axis regularity
4. 式 (4.8)–(4.11) の stress 定義を追加抽出
5. Sections 6–7 の pulse coordinate / phase / amplitude
6. Propositions 7.5–7.6 の stress realization

oscillatory realization は依存式が多いため、前段の数値部品を独立に検証してから進める。
