# Appendix B: Bρ ノルムと nonlinear fixed-point map 監査

## 対象

Proposition B.2 の係数空間 Bρ に対して、B.4 weight、B.7--B.10 型 convolution、nonlinear remainder `R1`,`R2`、実際の fixed-point map

`Phi -> Phi0 + (1+T)^-1 J2(R1)/(2 Lambda)`

`u -> u0 + J1(R2)/(2 Lambda)`

を監査する。

現在の full audit path では、fixed-point operator algebra の radial degree `alpha` と eta derivative order `beta` の双方を有限 cutoff に依存しない analytic majorant で評価する。また fixed multipliers に加え、center `u0=-Y Z*/(2L)` も complex neighborhood 上で閉じる。

## B.4 weight

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を使用する。finite jet 実装は回帰確認用として残すが、full certificate の定数決定には使わない。

## infinite alpha / beta convolution

B.7--B.10 の convolution mechanism を alpha と beta の両方に対して

`sum_{r=0}^n 1 / ((r+1)^2 (n-r+1)^2) <= C2/(n+1)^2`

`C2 = 4 pi^2 / 3`

で閉じる。`appendix_b_brho_infinite_alpha_beta_operator_audit(rho)` は radial cutoff を取らず、standalone の `D_X`, `I`, `J1`, `J2`, `d_eta` には有限 certified constant を与えない。実際の式に現れる composite operator だけを直接評価する。

`rho=1e-5` の代表値:

- `product_bound = 173.1717173937821`
- `mixed_AX_J2 = 8.8663919306e10`
- `pressure_J1_detaI = 7.6770415885e13`

これらは conservative reproduction majorant であり、原典の最適定数をそのまま転記したものではない。

## fixed multiplier の全 beta analytic bound

`chi`, `zeta*`, `1/L`, `H*`, `W*`, `g=phi*/C` は complex stadium

`Omega_R = { z : dist(z,[-1,1]) <= R }`

上の sup と Cauchy estimate で全 beta を評価する。ordered baseline `rho=1e-5` では fixed multiplier 側の `R=4e-5`、Cauchy factor は 1。

`g=phi*/C` は Lambda を先に選び、その後 C を十分大きく選んで `sup_Omega |g|<=1` とする。

## Pi0 と Z*/L の analytic center bound

`include/appendix_b_pressure_analytic.hpp` で Appendix A.4 pressure datum を complex neighborhood 上に持ち上げる。

`f(z)=1/(1+z^2)` の pole は `z=+/-i` で、`dist(z,[-1,1])<=R<1` なら

`|f(z)| <= (1-R)^-2`。

A.4 の eta 依存は次の構造を持つ。

- A.10 より前: 正の radial mass の定数倍 `f(eta)^2`
- A.10: `c(y) f(eta)^(2 theta(y))`, `0<=theta<=1`
- A.10 終了後: eta 依存が厳密に消える

したがって eta=0 の全正 radial mass を `sup |f|^2` 倍することで `Pi0/P_*^2` を外側 stadium 上で保守的に majorize できる。

baseline では

- pressure inner radius = `4e-5`
- pressure outer radius = `8e-5`
- `sup |f| <= 1.0001600192`
- `sup |Pi0|/P_*^2 <= 3.3156836263`
- inner stadium 上 `sup |d_eta Pi0|/P_*^2 <= 8.2892091e4`
- `||Pi0/P_*^2||_(Bρ) <= 3.3156836263`
- `||d_eta Pi0/P_*^2||_(Bρ) <= 4.1446045e4`

となる。

B.1 の

`Z* = -A(1-2 eta U*)U* - H* U*_eta - d Pi0_eta + 4 A eta Pi0`

を同じ inner stadium 上で評価し、`1/L` も掛けると

`||(Z*/L)/P_*^2||_(Bρ) <= 1.6579744445e5`

`log ||Z*/L||_(Bρ) <= 152.0185221081`

となる。よって `Y<=4.1` で

`u0 = -Y Z*/(2L)`

に対し

`||u0||_(Bρ) <= 2.1505976143e66`

を得る。これは旧 real-axis sampled center より約 `1.67e4` 倍保守的だが、有限 eta sampling への依存を完全に外している。

## (1+T)^-1

`T=(1/2)J2(chi .)` について

`K = 40 ||chi||_(Bρ)`

`||T^k|| <= K^k/[k!(k+1)!]`

を全 starting alpha に対して使う。`rho=1e-5` では

- `||chi||_(Bρ) <= 2.2960763861`
- `K = 91.8430554443`
- `||(1+T)^-1|| <= 1.9712281018e6`

となる。

## Proposition B.2 certificate

`AppendixBPropositionB2Certificate` は以下を同時に要求する。

- infinite alpha / beta operator certification
- complex-neighborhood fixed multiplier certification
- analytic center certification (`Pi0`, `Pi0_eta`, `Z*/L`, `u0`)
- inverse certification
- contraction `q<1`
- invariant ball `source + q r <= r`

analytic-center flag が false の場合は、他の条件が成立しても `proposition_b2_certificate=false` とする。

ordered baseline `rho=1e-5` では最初に certificate が閉じる decade は

`Lambda = 1e21 Lambda0`

で、そのとき

- `q = 0.64051349483873`
- source norm = `1.2625982006e6`
- ball radius = `3.6878383239e6`
- invariant LHS = `3.6247084138e6`
- invariant margin = `6.3129910031e4 > 0`
- `analytic_center_certified = 1`
- `infinite_alpha_beta = 1`
- `inverse_certified = 1`
- `contraction_certified = 1`
- `invariant_ball_certified = 1`
- `proposition_b2_certificate = 1`

となる。

これは最小 Lambda や原典の最適 constant を主張するものではなく、現在の conservative analytic majorant で Proposition B.2 の fixed-point 条件を明示的に閉じる一例である。

## 成果物

CI は以下を生成する。

- `results/reference/appendix_b_scale.csv`
- `results/reference/appendix_b_brho.csv`
- `results/reference/appendix_b_fixed_map.csv`
- `results/reference/appendix_b_prop_b2_certificate.csv`

certificate CSV には pressure radii、`Pi0/P_*^2`, `Pi0_eta/P_*^2`, `Z*/L`, `u0`, source、Lipschitz `q`、ball radius、invariant margin、全 certification flag を保存する。

## 現在の境界

今回、Proposition B.2 certificate から主要な有限切断・実軸中心依存

- `beta_max`
- `radial_order`
- sampled `max_eta |u0|`

を外した。

一方、原典証明の完全形式化とはまだ区別する。残る主な検証項目は次。

1. B.4 と B.7--B.10 の exact manuscript transcription、とくに numerical constants の逐語照合。
2. conservative slack constants の sharpening。
3. complex pole separation、polynomial sup、pressure mass majorant を interval arithmetic で独立 certification。
4. Proposition B.3 の positivity / endpoint shear inequality B.19 を、この full analytic profile certificate と接続する。
