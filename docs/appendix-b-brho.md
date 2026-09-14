# Appendix B: Bρ ノルムと nonlinear fixed-point map 監査

## 対象

Proposition B.2 の係数空間 Bρ に対して、B.4 weight、B.7--B.10 型 convolution、nonlinear remainder `R1`,`R2`、実際の fixed-point map

`Phi -> Phi0 + (1+T)^-1 J2(R1)/(2 Lambda)`

`u -> u0 + J1(R2)/(2 Lambda)`

を監査する。

現在の実装では、fixed-point map の operator algebra について radial degree `alpha` と eta derivative order `beta` の双方を有限 cutoff に依存しない analytic majorant へ置き換えた。

## B.4 weight

現在の転記では

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を使用する。

finite jet の実装は回帰確認用として残しているが、full audit path では fixed-point bound の定数決定には使わない。

## infinite alpha / beta convolution

B.7--B.10 の convolution mechanism を、alpha と beta の両方に対して同じ square-summable kernel で閉じる。

基本評価は

`sum_{r=0}^n 1 / ((r+1)^2 (n-r+1)^2)`

`<= C2 / (n+1)^2`

`C2 = 4 pi^2 / 3`

である。

B.4 では radial weight に `(alpha+1)^-2`、eta weight に `(beta+1)^-2` が入るため、product convolution は両方向を同じ型の bound で処理できる。

実装 `appendix_b_brho_infinite_alpha_beta_operator_audit(rho)` は radial cutoff を引数に取らない。旧 compatibility wrapper に `na=8` と `na=64` を与えても、product / mixed / pressure composite constant が完全に一致することを CI で回帰確認している。

full audit では standalone の

- `D_X`
- `I`
- `J1`
- `J2`
- `d_eta`

には有限 certified constant を割り当てず、意図的に `infinity` とする。

代わりに実際の方程式に現れる composite

- `Jnu(FG)`
- `Jnu(F D_X G)`
- `Jnu(F d_eta G)`
- `Jnu[(d_eta F)(D_X G)]`
- `A_X(F)` を含む同型 composite
- pressure の `J1 I`, `J1 d_eta I`, `J1 D_X I`

のみを直接評価する。

これは B.6--B.10 の cancellation / smoothing を operator の途中で壊さないためである。

### 現在の保守 majorant

`rho=1e-5` で現在使う代表値は

- `product_bound = 173.1717173937821`
- `mixed_AX_J2 = 8.8663919306e10`
- `pressure_J1_detaI = 7.6770415885e13`

である。

これらは `alpha=0..infinity`, `beta=0..infinity` を一括で含む reproduction majorant であり、有限 `radial_order` / `beta_max` の最大値探索から作っていない。

ただし、ここで使う数値 slack (`80`, `160`, eta-shift factor など) は原典の最適定数をそのまま転記したものではない。B.7--B.10 の構造に沿った保守的な再現 bound であり、primary manuscript の exact constant transcription は別の検証項目として残す。

## fixed multiplier の全 beta analytic bound

`chi`, `zeta*`, `1/L`, `H*`, `W*`, `g=phi*/C` は complex neighborhood 上の Cauchy estimate で全 beta を評価する。

`Omega_R = { z : dist(z,[-1,1]) <= R }`

上で `sup |f| <= M` なら

`|d_eta^beta f| <= M beta! R^(-beta)`

なので B.4 norm は

`||f||_(Bρ) <= M sup_beta (rho/R)^beta (beta+1)^2`

で閉じる。

ordered baseline `rho=1e-5` では `R=4e-5`、Cauchy factor は `1`。

`g=phi*/C` は原典どおり Lambda を先に選び、その後 C を十分大きく選ぶ。C 自体は巨大なので exponentiate せず `log_C_for_g` のみ保持する。

## (1+T)^-1

`T=(1/2)J2(chi .)` は radial degree を必ず 1 上げる。

`K = 40 ||chi||_(Bρ)`

として

`||T^k|| <= K^k / (k! (k+1)!)`

を全 starting alpha に対して使う。

したがって inverse の prefix length は radial cutoff ではなく、factorial series をどこで analytic tail に切り替えるかという単なる数値的 split である。

`rho=1e-5` では

- `||chi||_(Bρ) <= 2.2960763861`
- `K = 91.8430554443`
- `||(1+T)^-1|| <= 1.9712281018e6`

となる。

## all-alpha-beta fixed-point map

B.13 の `|U-U*|<=0.05` から得る基準は

`log Lambda0 = 145.2883536`

である。

`rho=1e-5`, `Phi_bound=2`, `u_bound=1.1 max|u0|` で、fixed multiplier と R1/R2 operator algebra の双方を infinite alpha / beta majorant にすると

- `Lambda = Lambda0`: contraction bound `2.0512241094e16`
- `Lambda = 1e16 Lambda0`: `2.0512241094`
- `Lambda = 1e18 Lambda0`: `0.02051224109` — contraction `<1`
- `Lambda = 1e20 Lambda0`: `2.0512241094e-4`

となる。

支配項は引き続き R1 の

`d zeta* u Phi`

channel である。

この Lambda factor は最小値でも原典の最適値でもない。現在の conservative analytic majorant で contraction を明示的に閉じる一例である。

## 成果物

CI は

- `results/reference/appendix_b_scale.csv`
- `results/reference/appendix_b_brho.csv`
- `results/reference/appendix_b_fixed_map.csv`

を生成する。

`appendix_b_fixed_map.csv` には

- `infinite_eta=1`
- `infinite_alpha=1`
- product / mixed / pressure composite constants
- complex-neighborhood multiplier bound
- `(1+T)^-1` factorial bound
- R1/R2 各 channel
- contraction certification

を保存する。

## 現在の境界

今回、fixed-point map の coefficient-space operator algebraから主要な有限切断

- `beta_max=6`
- `radial_order=18`

を外した。

したがって「有限 degree を大きくすれば閉じる」という数値監査から一段進み、現在の reproduction majorant の範囲では `alpha,beta` の全次数を同時に含む contraction audit になった。

一方、Proposition B.2 の原典証明を完全に形式化・再証明したとはまだ言わない。残る主な検証項目は次である。

1. B.4 と B.7--B.10 の exact manuscript transcription、とくに原典が用いる numerical constants を primary PDF と逐語的に照合する。
2. 現在の conservative slack constants を exact convolution algebra から sharpen し、必要 Lambda factor の過大評価を下げる。
3. complex pole separation / polynomial sup を interval arithmetic で独立 certification する。
4. Proposition B.2 の invariant-ball 条件を Lipschitz contraction と同じ full Bρ majorant で明示的にまとめる。

次の本筋は 4、すなわち **full Bρ invariant ball + contraction を一つの Proposition B.2 certificate に統合すること**である。
