# Appendix B: Bρ ノルムと nonlinear remainder / fixed-point map 監査

## 対象

Proposition B.2 の係数空間で用いられる (B.4) の重みと、(B.15) 後の nonlinear remainder `R1`,`R2`、さらに実際の fixed-point map

`Phi -> Phi0 + (1+T)^-1 J2(R1)/(2 Lambda)`

`u -> u0 + J1(R2)/(2 Lambda)`

を監査する。

ここでの計算は論文の無限次元 Banach-space 証明を置き換えない。現在は radial degree `alpha` について無限 factorial tail まで解析的に閉じ、さらに `chi`, `zeta*`, `1/L`, `H*`, `W*`, `g=phi*/C` など固定 eta multiplier は complex neighborhood 上の Cauchy estimate により **全 beta** を一括評価する。

ただし、積・B.6 mixed estimate 等の coefficient-space operator constant はまだ `beta_max=6` の有限和から評価している。そのため「multiplier の eta 方向」は全 beta 化したが、Proposition B.2 の Bρ operator algebra 全体が無限 beta で閉じたわけではない。

## B.4 weight

実装は論文の重み

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を直接使用する。

`include/appendix_b_brho.hpp` では、積の Leibniz 展開と `I`, `A_X`, `D_X`, `J1`, `J2`、B.6 型 composite operator の有限 `(alpha,beta)` constant を計算する。

## fixed multiplier の全 beta analytic bound

`include/appendix_b_analytic_multiplier.hpp` では、real-axis 上の有限 Taylor jet sampling の代わりに、`[-1,1]` の周囲の complex stadium

`Omega_R = { z : dist(z,[-1,1]) <= R }`

を使う。

scalar multiplier `f` が `Omega_R` で正則かつ

`sup_{Omega_R} |f| <= M`

なら、Cauchy estimate より

`|d_eta^beta f| <= M beta! R^(-beta)`

である。B.4 の `alpha=0` weight に割ると

`||f||_(Bρ) <= M * sup_{beta>=0} (rho/R)^beta (beta+1)^2`

となる。コードはこの最後の supremum を `scalar_cauchy_factor` として全 beta について評価する。

ordered baseline の `rho=1e-5` では `R=4e-5` とし、`rho/R=1/4` なので Cauchy factor は `1` である。

### H*, chi, zeta* の pole separation

`H*(z)` は

`H*(z) = j0 + (4.5-h) z - j0 z^2 - 4 z^3`

である。`chi` と `zeta*` の pole は `H*(z)=+/- i sigma*` にある。

数値 root finder は最寄り pole 距離の診断値を保存するが、正則性の保証そのものには使わない。実際の bound は `z=x+iy` に対する恒等式

`Im H*(x+iy) = y ( H*'(x) + 4 y^2 )`

を使う。

`G = sup_{x in [-1-R,1+R]} |H*'(x)|`

とすると `Omega_R` 上で

`|Im H*(z)| <= R (G + 4 R^2)`

である。そこで

`m = sigma* - R (G + 4 R^2)`

が正なら

`|H*(z) +/- i sigma*| >= m > 0`

となり、`Omega_R` に pole が存在しないことが直接分かる。

さらに

`chi = 1 - sigma*^2 / [(H*-i sigma*)(H*+i sigma*)]`

`H*/(H*^2+sigma*^2) = (1/2)[1/(H*-i sigma*) + 1/(H*+i sigma*)]`

を使うため、粗い `sup|H| / inf|H^2+sigma^2|` を使わず

`sup |chi| <= 1 + sigma*^2/m^2`

`sup |zeta*| <= sup|L|/m`

と評価できる。

### 1/L, H*, W*

`L(z)=1-2h z^2` については

`|L(z)| >= 1 - 2h(1+R)^2`

を使う。ordered baseline では `h=e^-80` なので 1/L はほぼ 1 のままである。

`H*`, `W*`, `eta`, `d`, `1-2eta U*`, R2 の linear-u coefficient, `4A eta` はすべて低次多項式なので、係数絶対値と `|z|<=1+R` から complex sup を直接上から抑え、その後 Cauchy factor を掛けて全 beta の Bρ norm とする。

### g=phi*/C

`phi*(z) = exp(Lambda integral_0^z zeta*(w) dw)`

である。Lambda が巨大なので `phi*` の高階微分を直接展開しない。論文の順序どおり Lambda を先に選び、その後 C を十分大きく取る。

`Omega_R` 内の 0 から z までの path length を `1+R` で抑えると

`log C >= Lambda (1+R) sup_{Omega_R}|zeta*|`

で

`sup_{Omega_R}|g| <= 1`

を保証できる。したがって g の全 beta derivative も同じ Cauchy estimate で閉じ、fixed-point budget には

`||g||_(Bρ) <= scalar_cauchy_factor`

を渡す。C 自体は桁が巨大なので exponentiate せず `log_C_for_g` のみ保存する。

## B.6 mixed estimate

`d_eta` と `D_X` を独立 operator として掛け合わせると非常に粗い bound になる。論文はこの組合せを

`J_nu[(d_eta F)(D_X G)]`

として一体で評価する。

実装では特に

- `J2[(d_eta A_X u)(D_X Phi)]`
- `J1[(d_eta A_X u)(D_X u)]`

に `mixed_AX_J2`, `mixed_AX_J1` を直接適用する。

現時点の重要な未完了点は、これら operator constants 自体がまだ `beta_max=6` の有限 convolution audit であること。fixed multiplier は全 beta 化済みだが、この operator algebra を B.7--B.10 の infinite-beta convolution inequality に置き換える必要がある。

## (1+T)^-1: finite prefix + analytic factorial tail

`T=(1/2)J2(chi .)` について単純な `||T||<1` は要求しない。`T` は radial degree を必ず1つ上げ、B.4 weight の比から

`step_alpha <= K / ((alpha+1)(alpha+2))`

`K = 40 ||chi||_(Bρ)`

を使う。start degree 0 では

`||T^k|| <= K^k / (k! (k+1)!)`

となる。

`include/appendix_b_full_inverse.hpp` は有限 k prefix と analytic factorial tail を組み合わせる。fixed multiplier を全 beta Cauchy bound に置き換えた後は `||chi||_(Bρ)` が保守的に大きくなるため、`include/appendix_b_fixed_map_budget.hpp` は finite prefix の終了 k を固定 18 にせず、最初の omitted factorial ratio が 1/2 未満になるまで自動延長する。

ordered baseline `rho=1e-5` の all-beta multiplier audit では

- `analytic radius R = 4e-5`
- 最寄り `H*=+/-i sigma*` root の診断距離: 約 `5.556e-4`
- Cauchy factor: `1`
- `||chi||_(Bρ) <= 2.2960763861`
- `||zeta*||_(Bρ) <= 455.38140254`
- `K = 91.84305544`
- full radial inverse bound: `10587.80005495`
  - finite prefix: `10501.93259851`
  - factorial tail: `85.86745644`

となった。

## analytic eta multiplier を入れた fixed-point map

B.13 の `|U-U*|<=0.05` から得る基準値は

`log Lambda0 = 145.2883536`

である。

`rho=1e-5`, `Phi_bound=2`, `u_bound=1.1 max|u0|` で、fixed multiplier と g を all-beta complex-neighborhood bound に置き換えた結果、fixed-point Lipschitz bound は

- `Lambda = Lambda0`: `1.523708526e8`
- `Lambda = 1e6 Lambda0`: `152.3708526`
- `Lambda = 1e7 Lambda0`: `15.23708526`
- `Lambda = 1e8 Lambda0`: `1.523708526`
- `Lambda = 1e9 Lambda0`: `0.1523708526` — contraction `<1`

となる。

支配項は引き続き R1 の

`d zeta* u Phi`

を通る channel である。all-beta bound は finite-beta sampled multiplier よりかなり保守的なので、必要な Lambda decade は `1e7 Lambda0` から `1e9 Lambda0` へ上がった。この値は論文の最適定数や最小 Lambda を意味しない。

`log_C_for_g` は Lambda に比例し、基準 Lambda0 でも非常に大きい。これは C を数値的に実体化する必要があるという意味ではなく、原典の「Lambda を選んだ後 C>=C0(Lambda) を選ぶ」という存在証明の順序を log-space で監査している。

## 成果物

CI は次を生成する。

- `results/reference/appendix_b_scale.csv`: B.13 から得る ordered Lambda scale
- `results/reference/appendix_b_brho.csv`: rho sweep、operator constants、all-beta analytic multiplier norms、raw R1/R2 budget
- `results/reference/appendix_b_fixed_map.csv`: complex radius / pole-separation diagnostics / Cauchy factor / `log_C_for_g` / full radial inverse / fixed-map 各項 / contraction certification

## 未完了境界

今回、requested multiplier

- `chi`
- `zeta*`
- `1/L`
- `H*`
- `W*`
- `g=phi*/C`

の eta 方向は complex neighborhood 上の all-beta Cauchy bound へ置き換えた。

Proposition B.2 の full Bρ analytic contraction と呼ぶには、次が残る。

1. `appendix_b_brho.hpp` の product / B.6 mixed / Jnu composite constant を `beta_max=6` の有限列挙から、B.7--B.10 に基づく infinite-beta convolution bound へ置き換える。
2. pressure channel の `J1 d_eta I` / `J1 D_X I` を専用 composite bound に置き換え、現在の保守的 standalone bound を削る。
3. rho と radial finite-prefix split を変えたときの analytic bound / 必要 Lambda factor の安定性を確認する。
4. 数値 floating-point root 距離は cross-check に留めているが、必要なら interval arithmetic による独立検証を追加する。

radial alpha tail と固定 multiplier の eta derivative tail は切断依存を除いた。次の主要な切断依存は B.7--B.10 operator convolution の beta 方向である。
