# Appendix B: Bρ ノルムと nonlinear remainder / fixed-point map 監査

## 対象

Proposition B.2 の係数空間で用いられる (B.4) の重みと、(B.15) 後の nonlinear remainder `R1`,`R2`、さらに実際の fixed-point map

`Phi -> Phi0 + (1+T)^-1 J2(R1)/(2 Lambda)`

`u -> u0 + J1(R2)/(2 Lambda)`

を監査する。

ここでの計算は論文の無限次元 Banach-space 証明を置き換えない。現在は radial degree `alpha` について無限 tail まで解析的に閉じたが、eta derivative は `beta<=6` の Taylor jet と実軸上でサンプリングした multiplier norm を使用している。したがって、以下でいう full inverse は **radial alpha 方向の full tail** を意味し、complex neighborhood 上の全 beta を含む Proposition B.2 全体の証明を意味しない。

## B.4 weight

実装は論文の重み

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を直接使用する。

有限 beta 監査では

`max |d_eta^beta F_alpha| / a(alpha,beta)`

を評価する。`include/appendix_b_brho.hpp` では、積の Leibniz 展開を有限和で全列挙し、`I`, `A_X`, `D_X`, `J1`, `J2` と composite operator の bound を計算する。

## B.6 mixed estimate

`d_eta` と `D_X` を独立 operator として掛け合わせると非常に粗い bound になる。論文はこの組合せを

`J_nu[(d_eta F)(D_X G)]`

として一体で評価する。

実装ではさらに、derivative を省いた版も同じ有限和から評価し、特に fixed-point map では

- `J2[(d_eta A_X u)(D_X Phi)]`
- `J1[(d_eta A_X u)(D_X u)]`

に `mixed_AX_J2`, `mixed_AX_J1` を直接適用する。standalone `d_eta * D_X` bound はこの2項には使わない。

## fixed-point map の term-wise budget

`include/appendix_b_fixed_map_budget.hpp` では、`R1/R2` の raw norm を先に1個の M に潰さず、`J2/(2 Lambda)` または `J1/(2 Lambda)` を含めた各項の Lipschitz contribution を保存する。

`W` も

`W = W* + Lambda^-1[-2D eta A_X(u) - d d_eta A_X(u)]`

まで分解するため、B.6 の mixed structure を保持できる。

## (1+T)^-1: finite prefix + analytic factorial tail

`T=(1/2)J2(chi .)` について単純な `||T||<1` は要求しない。ordered baseline、`rho=1e-5` では crude one-step bound は

`||T|| <= 19.999964`

であり、幾何級数 `1/(1-||T||)` は使用不能である。

一方、`T` は radial degree を必ず1つ上げる。B.4 weight の比から各 radial degree `alpha` について

`step_alpha <= K / ((alpha+1)(alpha+2))`

という envelope を使う。現在の実装では

`K = 40 ||chi||_(Bρ)`

で、`rho=1e-5` の ordered baseline では

`K = 39.99992839296`

となる。したがって start degree 0 からは

`||T^k|| <= K^k / (k! (k+1)!)`

となり、`||T||>1` でも Neumann series の radial tail は factorial に収束する。

実装 `include/appendix_b_full_inverse.hpp` は次の hybrid bound を使う。

1. `k=0..18` は B.4 の exact finite-beta step constant を積み上げる。
2. start radial degree `s=0..64` は直接走査する。
3. 未走査の `s>=65` は `K/((s+j+1)(s+j+2))` envelope で全て包含する。
4. `k>=19` は `K^k/(k!(k+1)!)` の analytic tail を使う。
5. tail の連続項比が単調減少することを使い、最初の omitted term から残り全部を geometric majorant で閉じる。

`rho=1e-5`, `beta_max=6` の結果は

- finite prefix bound: `352.87245276559923`
- analytic factorial tail: `1.0265295066004976e-5`
- full radial inverse bound: `352.87246303089432`
- envelope regression ratio: `1.0`

となった。

B.4 weight の `20^-alpha` は direct double では非常に大きい alpha で underflow するため、envelope の実装回帰チェックは `alpha<=128` に制限する。これは無限 alpha の保証を有限 scan に置き換えるものではない。無限 alpha は上記 algebraic envelope で抑え、有限 scan はコード実装の回帰検査としてのみ使う。

## rho の選択

ordered baseline では separation audit が `sigma*=0.0025` を与える。`zeta*=-L H*/(H*^2+sigma*^2)` の最も近い複素 pole の距離は概ね `sigma*/|H*'(eta0)|` で、今回の値では約 `5.6e-4` である。

したがって `rho=0.01` は analytic neighborhood として大きすぎる。診断では `rho=1e-5 ... 2e-4` を中心に sweep する。

## full radial inverse を入れた fixed-point map の結果

`rho=1e-5`, `beta_max=6`, `Phi_bound=2`, `u_bound=1.1 max|u0|`, `||g||<=1` の監査で、B.13 の `|U-U*|<=0.05` だけから得る基準値は

`log Lambda0 = 145.2883536`

である。

full radial inverse bound `352.8724630` を fixed-point map に流すと、Lipschitz bound は

- `Lambda = Lambda0`: `1.115119e6`
- `Lambda = 1e5 Lambda0`: `11.15119`
- `Lambda = 1e6 Lambda0`: `1.115119` — まだ contraction ではない
- `Lambda = 1e7 Lambda0`: `0.1115119` — contraction `<1`

となる。

支配項は引き続き `R1` の `d zeta* u Phi` を通る

`(1+T)^-1 J2[d zeta* u Phi / L] / (2 Lambda)`

であり、B.6 mixed derivative term は支配項ではない。

旧 finite-degree inverse (`119.14753`) では `1e6 Lambda0` で contraction が閉じていたが、無限 radial tail を含む保守的 boundへ置き換えると必要 decade は `1e7 Lambda0` に1桁上がった。これは期待される方向の変化である。

以上により、**finite beta の B.4 audit の範囲では、radial degree を無限まで含めた `(1+T)^-1` bound と実際の fixed-point map contraction を接続できた**。

## 成果物

CI は次を生成する。

- `results/reference/appendix_b_scale.csv`: B.13 から得る ordered `Lambda` scale
- `results/reference/appendix_b_brho.csv`: rho sweep、B.4/B.6 operator constants、multiplier norms、raw R1/R2 budget
- `results/reference/appendix_b_fixed_map.csv`: full radial inverse の finite prefix / analytic tail / K / envelope ratio、fixed-map 各項、`M_Phi_map`, `M_u_map`, contraction certification

## 未完了境界

Proposition B.2 の full analytic contraction と呼ぶには、まだ次が必要。

1. `chi`, `zeta*`, `1/L`, `H*`, `W*` 等の multiplier norm を、有限 `beta<=6` の real-axis sample ではなく complex neighborhood 上の全 beta の解析上界へ置き換える。
2. `g=phi*/C` の complex-neighborhood bound を C の選択とともに閉じる。
3. pressure channel の `J1 d_eta I` / `J1 D_X I` を専用 composite bound に置き換え、現在の保守的 standalone bound を削る。
4. beta order と rho を変えたときの fixed-map bound / 必要 Lambda factor の安定性を確認する。

radial alpha tail 自体は今回の finite-prefix + analytic-factorial-tail で切断依存を除いた。
