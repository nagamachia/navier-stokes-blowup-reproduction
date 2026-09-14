# Appendix B: Bρ ノルムと nonlinear remainder / fixed-point map 監査

## 対象

Proposition B.2 の係数空間で用いられる B.4 の重みと、B.15 後の nonlinear remainder `R1`,`R2`、さらに実際の fixed-point map

`Phi -> Phi0 + (1+T)^-1 J2(R1)/(2 Lambda)`

`u -> u0 + J1(R2)/(2 Lambda)`

を監査する。

現在の実装では

- radial degree `alpha` の `(1+T)^-1` tail は factorial majorant で無限まで閉じる。
- `chi`, `zeta*`, `1/L`, `H*`, `W*`, `g=phi*/C` 等の固定 eta multiplier は complex neighborhood 上の Cauchy estimate により全 beta を一括評価する。
- product / B.6 mixed / Jnu composite / pressure composite は B.7--B.10 型の beta-independent convolution majorant に置き換え、fixed-point map の実経路から `beta_max=6` 依存を除いた。

ただし、以下の all-beta 数値定数は原典の構造を再現するための保守的 majorant であり、B.7--B.10 の最適定数または原典記載の定数をそのまま転記したものとは扱わない。primary manuscript の exact typesetting / constant の再確認は別途必要である。

## B.4 weight

現在の転記では

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を使用する。

## fixed multiplier の全 beta analytic bound

`include/appendix_b_analytic_multiplier.hpp` では、`[-1,1]` の周囲の complex stadium

`Omega_R = { z : dist(z,[-1,1]) <= R }`

を使う。scalar multiplier `f` が `Omega_R` で正則かつ

`sup_{Omega_R} |f| <= M`

なら Cauchy estimate より

`|d_eta^beta f| <= M beta! R^(-beta)`

である。B.4 の `alpha=0` weight に割ると

`||f||_(Bρ) <= M * sup_{beta>=0} (rho/R)^beta (beta+1)^2`

となる。

ordered baseline の `rho=1e-5` では `R=4e-5` とし、Cauchy factor は `1` である。

### H*, chi, zeta* の pole separation

`H*(z) = j0 + (4.5-h) z - j0 z^2 - 4 z^3`

で、`chi` と `zeta*` の pole は `H*(z)=+/- i sigma*` にある。

数値 root finder は最寄り pole 距離の cross-check にのみ使い、正則性の保証には

`Im H*(x+iy) = y (H*'(x) + 4 y^2)`

を使う。`G=sup |H*'|` として

`m = sigma* - R (G + 4R^2) > 0`

なら

`|H*(z) +/- i sigma*| >= m`

である。

さらに

`chi = 1 - sigma*^2 / [(H*-i sigma*)(H*+i sigma*)]`

`H*/(H*^2+sigma*^2) = (1/2)[1/(H*-i sigma*) + 1/(H*+i sigma*)]`

を使い、`chi` と `zeta*` を直接評価する。

### g=phi*/C

`phi*(z)=exp(Lambda integral_0^z zeta*(w)dw)` に対して、原典の順序どおり Lambda を選んだ後 C を選ぶ。

`log C >= Lambda (1+R) sup_{Omega_R}|zeta*|`

とすれば `sup_{Omega_R}|g|<=1` であり、g の全 beta derivative も Cauchy estimate で閉じる。C は exponentiate せず `log_C_for_g` のみ保存する。

## B.7--B.10 型 infinite-beta convolution majorant

`include/appendix_b_brho.hpp` の `appendix_b_brho_infinite_beta_operator_audit` が all-beta proof path を担当する。

beta 重みの基本畳み込みについて

`sum_{r=0}^beta 1 / ((r+1)^2 (beta-r+1)^2)`

を beta に依存しない形で

`<= C2 / (beta+1)^2`

と抑え、実装では明示的な保守定数

`C2 = 4 pi^2 / 3`

を使う。product majorant は

`C_alg = C2^2`

とし、radial index と eta derivative shift は別の explicit factor で吸収する。

現在 `radial_order=18` の audit では

- `product_bound = 173.1717173937821`
- `mixed_AX_J2 = 3.9898763687527388e14` at `rho=1e-5`
- `pressure_J1_detaI = 3.5465567722246573e17` at `rho=1e-5`

となる。

これはかなり保守的だが、beta=0..infinity を一括して含む majorant である。

### standalone operator を使わない理由

all-beta path では standalone の

- `I`
- `J1`
- `J2`
- `d_eta`

には有限 certified constant を割り当てず、意図的に `infinity` とする。

代わりに、実際の方程式に現れる composite

- `Jnu(FG)`
- `Jnu(F D_X G)`
- `Jnu(F d_eta G)`
- `Jnu[(d_eta F)(D_X G)]`
- `A_X(F)` を含む同型の composite
- pressure の `J1 I`, `J1 d_eta I`, `J1 D_X I`

だけを直接 all-beta で抑える。

これは、`d_eta` や `D_X` を独立 operator norm に分離して巨大化させず、B.6--B.10 の構造を保持するためである。また old finite-beta path が誤って使われた場合には `infinity` が伝播して CI で露出する。

## (1+T)^-1 の all-beta inverse

`T=(1/2)J2(chi .)` について

`K = 40 ||chi||_(Bρ)`

とし、全 beta に対する factorial envelope

`||T^k|| <= K^k / (k! (k+1)!)`

を直接使う。

`appendix_b_full_inverse_all_beta_audit` は有限 eta jet を一切参照せず、tail ratio が 1/2 未満になるまで k-prefix を積み、その後を geometric majorant で閉じる。

`rho=1e-5` では

- `||chi||_(Bρ) <= 2.2960763861`
- `K = 91.8430554443`
- all-beta inverse bound `||(1+T)^-1|| <= 1.9655515236e6`
- prefix `1.9654656562e6`
- tail `85.86745644`

となる。

## all-beta fixed-point map の結果

B.13 の `|U-U*|<=0.05` から得る基準値は

`log Lambda0 = 145.2883536`

である。

`rho=1e-5`, `Phi_bound=2`, `u_bound=1.1 max|u0|` で、fixed multiplier と operator algebra の双方を all-beta majorant にした fixed-point Lipschitz bound は

- `Lambda = Lambda0`: `1.0226585827e18`
- `Lambda = 1e18 Lambda0`: `1.0226585827` — わずかに 1 より大きい
- `Lambda = 1e20 Lambda0`: `0.01022658583` — contraction `<1`
- `Lambda = 1e24 Lambda0`: `1.0226586e-6`

となる。

支配項は引き続き R1 の

`d zeta* u Phi`

を通る channel であり、pressure channel はこの ordered baseline では支配的ではない。

必要 Lambda factor が finite-beta audit より大きくなったのは、今回の all-beta operator constants が proof-oriented な保守 majorantだからである。これは最小 Lambda や原典の最適定数を意味しない。

rho sweep では、現在の保守 bound の範囲で

- `rho=1e-5`: `1e20 Lambda0` で閉じる。
- `rho=2e-5`: `1e20 Lambda0` で閉じる。
- `rho=5e-5`: 約 `1e28 Lambda0` で閉じる。
- `rho=1e-4`: `1e40 Lambda0` で閉じる。
- `rho=2e-4`: 現在の `1e40` sweep では閉じない。

したがって、今回の conservative analytic majorant では小さい rho の方が有利である。

## CI / 成果物

CI は次を生成する。

- `results/reference/appendix_b_scale.csv`
- `results/reference/appendix_b_brho.csv`
- `results/reference/appendix_b_fixed_map.csv`

`appendix_b_fixed_map.csv` には

- `infinite_eta=1`
- all-beta product / mixed / pressure composite constants
- complex radius / pole-separation diagnostics
- Cauchy factor / `log_C_for_g`
- all-beta `(1+T)^-1` bound
- fixed-map 各項
- contraction certification

を保存する。

## 現在の未完了境界

今回、eta derivative beta 方向については

- fixed multiplier
- product algebra
- B.6 mixed term
- Jnu composites
- pressure composites
- `(1+T)^-1`

から `beta_max=6` の切断依存を除いた。

ただし Proposition B.2 の full Banach-space 証明再現と呼ぶには、まだ次が必要である。

1. B.4 および B.7--B.10 の exact manuscript transcription / constant を primary PDF と照合する。現在の all-beta constants は保守的 reproduction majorant であり、原典の最適定数とは主張しない。
2. nonlinear product / composite operator の radial alpha 方向はまだ `radial_order=18` 由来の polynomial majorantを使う。`(1+T)^-1` の radial tail は無限まで閉じたが、R1/R2 operator algebraの alpha 方向は次の主要な有限切断である。
3. 必要なら interval arithmetic を導入し、complex pole separation / polynomial sup を floating-point cross-checkから独立に検証する。
4. constants を sharpen して必要 Lambda factor の保守性を下げる。

次の本筋は、B.7--B.10 の radial convolution についても `alpha=0..infinity` の explicit bound を作り、`radial_order=18` 依存を除くことである。
