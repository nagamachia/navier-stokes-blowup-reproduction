# Appendix B: Bρ ノルムと nonlinear remainder / fixed-point map 監査

## 対象

Proposition B.2 の係数空間で用いられる (B.4) の重みと、(B.15) 後の nonlinear remainder `R1`,`R2`、さらに実際の fixed-point map

`Phi -> Phi0 + (1+T)^-1 J2(R1)/(2 Lambda)`

`u -> u0 + J1(R2)/(2 Lambda)`

を有限 `(alpha,beta)` 切断で監査する。

ここでの計算は論文の無限次元 Banach-space 証明を置き換えない。特に、実軸上でサンプリングした multiplier norm と有限次数の eta Taylor jet を使う部分は数値監査であり、complex neighborhood 全体の厳密上界ではない。

## B.4 weight

実装は論文の重み

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を直接使用する。

有限切断ノルムは

`max |d_eta^beta F_alpha| / a(alpha,beta)`

で評価する。`include/appendix_b_brho.hpp` では、積の Leibniz 展開を有限和で全列挙し、product operator constant を計算する。また `I`, `A_X`, `D_X`, `J1`, `J2` の coefficient action から有限切断 operator bound を直接計算する。

## B.6 mixed estimate

`d_eta` と `D_X` を独立 operator として掛け合わせると非常に粗い bound になる。論文はこの組合せを (B.6)

`J_nu[(d_eta F)(D_X G)]`

として一体で評価する。

実装ではさらに、論文の「derivatives omitted」に対応して

- `J_nu[F G]`
- `J_nu[F D_X G]`
- `J_nu[F d_eta G]`
- `J_nu[(d_eta F)(D_X G)]`

を同じ有限和から直接評価する。`F` を `A_X(F)` に置き換える版も実装した。

これにより、fixed-point map では特に

- `J2[(d_eta A_X u)(D_X Phi)]`
- `J1[(d_eta A_X u)(D_X u)]`

に `mixed_AX_J2`, `mixed_AX_J1` を直接適用し、standalone `d_eta * D_X` bound は使用しない。

## R1 / R2 の term-wise M

raw remainder については `include/appendix_b_remainder_budget.hpp` で、原典の表示式に対応して各 contribution と `M_R1`, `M_R2` を保存する。

一方、実際の固定点写像については `include/appendix_b_fixed_map_budget.hpp` で、`J2/(2 Lambda)` または `J1/(2 Lambda)` を各項に先に適用した後の Lipschitz contribution を保存する。

`W` は一括評価せず、

`W = W* + Lambda^-1[-2D eta A_X(u) - d d_eta A_X(u)]`

まで分解する。そのため B.6 の mixed structure を失わない。

## (1+T)^-1

`T=(1/2)J2(chi .)` について単純な `||T||<1` は要求しない。

ordered baseline、`rho=1e-5`、有限切断 `(alpha_max,beta_max)=(18,6)` では crude one-step bound は

`||T|| <= 19.99996`

となり、幾何級数 `1/(1-||T||)` は使用不能である。

しかし `T` は radial degree を必ず1つ上げる。そこで各 radial degree `alpha` ごとに B.4 weight と eta Leibniz 展開を使って1-step ratio を計算し、出力 degree ごとに有限 Neumann sum

`1 + T + T^2 + ... + T^alpha`

の絶対値 bound を積み上げる。この有限切断 degree-raising bound では

`||(1+T)^-1|| <= 119.14753`

となった。

これは論文の factorial degree-raising mechanism の有限切断版である。無限次元 tail 全体の厳密 bound はまだ別途必要。

## rho の選択

ordered baseline では separation audit が `sigma*=0.0025` を与える。`zeta*=-L H*/(H*^2+sigma*^2)` の最も近い複素 pole の距離は概ね `sigma*/|H*'(eta0)|` で、今回の値では約 `5.6e-4` である。

したがって `rho=0.01` は analytic neighborhood として大きすぎる。原典と同じ順序、すなわち `sigma*` を固定した後でさらに小さい `rho` を選ぶ必要がある。診断では `rho=1e-5 ... 2e-4` を中心に sweep する。

## fixed-point map の結果

`rho=1e-5`, `(alpha_max,beta_max)=(18,6)`、`Phi_bound=2`, `u_bound=1.1 max|u0|`, `||g||<=1` の有限切断監査では、B.13 の `|U-U*|<=0.05` だけから得る基準値

`log Lambda0 = 145.28835`

に対し、実際の fixed-point map の Lipschitz bound は

`max(M_Phi_map,M_u_map) = 3.7652e5`

となった。

支配項は `R1` の `d zeta* u Phi` を通る

`(1+T)^-1 J2[d zeta* u Phi / L] / (2 Lambda)`

である。

Lambda は Appendix B の先行選択後にさらに十分大きく取れるため、同じ `rho`, `sigma*`, `P*` 等を固定したまま Lambda を増やすと contraction bound は低下する。今回の decade sweep では

- `Lambda = Lambda0`: `3.7652e5`
- `Lambda = 1e5 Lambda0`: `3.7652`
- `Lambda = 1e6 Lambda0`: `0.37652` — contraction `<1`
- `Lambda = 1e7 Lambda0`: `0.037652`

となった。

したがってこの有限切断モデルでは、B.6 を各 mixed term に直接適用し、degree-raising `(1+T)^-1` bound を使うことで、**十分大きい Lambda に対して実際の fixed-point map が contraction になることを数値的に再現できた**。

これは Proposition B.2 の無限次元証明そのものではない。特に必要な Lambda factor `1e6` は今回の保守的 finite-truncation / sampled-multiplier bound に依存し、論文中の最適な定数を意味しない。

## 成果物

CI は次を生成する。

- `results/reference/appendix_b_scale.csv`: B.13 から得る ordered `Lambda` scale
- `results/reference/appendix_b_brho.csv`: rho sweep、B.4/B.6 operator constants、multiplier norms、raw R1/R2 budget
- `results/reference/appendix_b_fixed_map.csv`: `rho`, Lambda factor、`chi/zeta` norm、`T` bound、degree-raising inverse bound、fixed-map 各項、`M_Phi_map`, `M_u_map`, contraction certification

## 未完了境界

次を満たすまでは Proposition B.2 の full analytic contraction を完了扱いにしない。

1. finite `(alpha,beta)` degree-raising bound を論文の infinite-dimensional factorial tail estimate へ接続する。
2. `g=phi*/C` の complex-neighborhood bound を C の選択とともに閉じる。
3. pressure channel の `J1 d_eta I` / `J1 D_X I` も専用 composite bound に置き換え、現在の保守的 standalone bound を削る。
4. radial/eta truncation を増やしたときの fixed-map bound と必要 Lambda factor の安定性を確認する。
