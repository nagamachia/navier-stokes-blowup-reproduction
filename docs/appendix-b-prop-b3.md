# Appendix B: Proposition B.3 継承 certificate

## 目的

Proposition B.2 で得た同じ full `B_rho` invariant ball から、Proposition B.3 の

- `Phi > 0` on `0 <= Y <= 4.1`
- endpoint `Y=4` の shear inequality (B.19)

を別の有限 radial truncation や独立な correction sweep を導入せずに継承する。

ここでいう certificate は現在の reproduction majorant に対する監査であり、原典の形式証明そのものではない。

## B_rho norm から pointwise bound への移行

B.4 の `beta=0` weight は

`a(alpha,0) = 20^(-alpha)/(alpha+1)^2`

なので、correction `delta Phi = Phi-Phi0` が

`||delta Phi||_(B_rho) <= r`

を満たせば、`0 <= Y < 20` について

`|delta Phi(Y,eta)| <= r S0(Y)`

`|D_X delta Phi(Y,eta)| <= r S1(Y)`

ただし

`S0(Y) = sum_{alpha>=0} (Y/20)^alpha/(alpha+1)^2`

`S1(Y) = sum_{alpha>=0} alpha (Y/20)^alpha/(alpha+1)^2`

である。

この和は `alpha=0..infinity` を直接含み、radial cutoff は再導入していない。

基準値は

- `S0(4.1) = 1.0565403893938534`
- `S0(4.0) = 1.0550188771985236`
- `S1(4.0) = 0.060698879372524904`

である。

## Phi > 0 の継承

B.11 の基準解について

`Phi0(Y,eta) = f0(Y chi(eta))`

かつ `0 <= Y chi <= 4.1` で

`f0 > 0.265`

を使う。

したがって B.2 ball radius `r` から

`Phi >= 0.265 - r S0(4.1)`

を得る。

B.2 certificate 自体は `Lambda = 1e21 Lambda0` から閉じるが、この時点では ball radius が大きすぎるため positivity はこの単純な norm inheritance では結論できない。

同じ certificate family で Lambda を増やすと `r` は縮み、最初に B.3 まで同時に閉じる decade は

`Lambda approximately 1e28 Lambda0`

である。

この点では

- `r = 0.13257281955643541`
- `r S0(4.1) = 0.14006853839719732`
- `Phi >= 0.12493146160280269 > 0`

となる。

## endpoint shear: chi > .99 branch

既存の B.19 baseline audit は `chi>.99` 上で unperturbed endpoint quantity

`p1(Phi0) = -2 D_X Phi0 / Phi0`

に対して

`p1(Phi0) >= 3.326818903256771`

を与える。

B.2 correction に対して

`e0 = r S0(4)`

`e1 = r S1(4)`

と置くと、`Phi0>=0.265` を用いた保守的な quotient bound

`p1(Phi) >= [p1(Phi0)*0.265 - 2 e1]/[0.265 + e0]`

を使える。

`Lambda approximately 1e28 Lambda0` では

- `e0 = 0.13986682723547295`
- `e1 = 0.008047021582331584`
- `p1(Phi) >= 2.1377719980377492`

となる。

よって azimuthal branch だけでも `>2` に明示的な余裕がある。

## endpoint shear: axial branch

`chi<=.99` 側は既存 B.19 baseline audit の `|Z*|` separation を継承する。

`Phi` perturbation の relative error

`epsilon = e0/0.265`

を吸収するため、既存の logarithmic `C` threshold を

`log C_required = log C_baseline - log(1-epsilon)`

まで保守的に増やす。

一方 `g=phi*/C` の complex-neighborhood normalization で選ぶ `C` は同じ Lambda とともに増大するため、上記 threshold を十分上回る。

最初の B.3 certificate 点では

- `log C chosen = 5.7026137929827761e93`
- `log C required = 8.3398243375995614e90`
- margin = `5.6942739686451762e93 > 0`

であり axial branch も閉じる。

## B.19 certificate

二つの endpoint branch を合わせ、現在の reproduction audit では

`c_ex = 0.1377719980377492`

を取れる。

したがって原典の endpoint shear 条件に対応する正の余裕を保持したまま、B.2 invariant ball から B.3 へ連続的に継承できる。

`results/reference/appendix_b_prop_b3_certificate.csv` には各 Lambda decade ごとに

- B.2 certificate
- ball radius
- `S0`,`S1`
- Phi lower bound
- endpoint correction errors
- unperturbed / perturbed shear lower bound
- `C` threshold / margin
- positivity / branch / B.19 / B.3 flags

を記録する。

## 現在の境界

B.2 から correction を B.3 へ運ぶ部分は full `B_rho` ball と infinite-alpha point-evaluation sum で処理されている。

ただし B.19 の **baseline endpoint branch audit** にはまだ real-eta sampling が残る。`appendix_b_endpoint_asymptotic_audit` は uniform eta samples に加えて狭い axial branch を落とさないため `eta0` を明示的に評価している。

したがって現在の正確な結論は、

**既存の B.19 baseline audit を anchor として、full B_rho Proposition B.2 certificate の perturbation を Proposition B.3 の positivity と endpoint inequality まで切断なしに継承した**

というものである。

次の改善点は、`chi>.99` / `|Z*|>delta*` の endpoint branch split 自体を interval / complex-neighborhood bound で全 eta certification し、B.19 baseline から sampling 依存も除くことである。
