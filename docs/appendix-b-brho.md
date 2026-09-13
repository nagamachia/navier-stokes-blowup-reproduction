# Appendix B: Bρ ノルムと nonlinear remainder 監査

## 対象

Proposition B.2 の係数空間で用いられる (B.4) の重みと、(B.15) 後の nonlinear remainder `R1`,`R2` を有限 `(alpha,beta)` 切断で監査する。

ここでの計算は論文の無限次元 Banach-space 証明を置き換えない。特に、実軸上でサンプリングした multiplier norm と有限次数の eta Taylor jet を使う部分は数値監査であり、complex neighborhood 全体の厳密上界ではない。

## B.4 weight

実装は論文の重み

`a(alpha,beta) = 20^(-alpha) rho^(-beta) beta! binom(alpha+beta,beta) / ((alpha+1)^2 (beta+1)^2)`

を直接使用する。

有限切断ノルムは

`max |d_eta^beta F_alpha| / a(alpha,beta)`

で評価する。`include/appendix_b_brho.hpp` では、積の Leibniz 展開を有限和で全列挙し、product operator constant を計算する。また `I`, `A_X`, `D_X`, `J1`, `J2` の coefficient action から有限切断 operator bound を直接計算する。

## B.6 mixed estimate

`d_eta` と `D_X` を独立 operator として掛け合わせると非常に粗い bound になる。実際、`rho=0.01`, radial order 18, eta order 6 の試行では standalone `d_eta` bound が約 `1.76e3` となり、raw contraction proxy は閉じなかった。

論文はこの組合せを (B.6)

`J_nu[(d_eta F)(D_X G)]`

として一体で評価する。このため `appendix_b_brho_mixed_Jnu_bound` を追加し、`F` および `A_X(F)` の両方について有限和から mixed constant を直接計算する。

## R1 / R2 の term-wise M

`include/appendix_b_remainder_budget.hpp` では、Lipschitz budget をブラックボックスの1定数にせず、原典の表示式に対応して分解する。

R1:

- `[W + h(1-2 eta U) + d u zeta*] Phi`
- `W D_X Phi`
- `Hc d_eta Phi`

R2:

- `[A(1-4 eta U*) + d U*_eta] u`
- `-2 A eta Lambda^-1 u^2`
- `W D_X u`
- `H* d_eta u`
- `d Lambda^-1 u d_eta u`
- `-4 A eta p`
- `d d_eta p`
- `-2 eta D_X p`

各 contribution と `M_R1`, `M_R2`, `M=max(M_R1,M_R2)` を CSV に保存する。

## rho の選択

ordered baseline では separation audit が `sigma*=0.0025` を与える。`zeta*=-L H*/(H*^2+sigma*^2)` の最も近い複素 pole の距離は概ね `sigma*/|H*'(eta0)|` で、今回の値では約 `5.6e-4` である。

したがって `rho=0.01` は analytic neighborhood として大きすぎる。原典と同じ順序、すなわち `sigma*` を固定した後でさらに小さい `rho` を選ぶ必要がある。診断では `rho=1e-5 ... 5e-4` を sweep する。

## 最初の非閉包から得た情報

不適切な `rho=0.01` で raw bound を計算したとき、CI は意図的に非閉包を検出した。

- `log Lambda = 145.288`
- `log M = 174.719`
- raw `log contraction proxy = 26.3022`
- `||zeta*||_(Bρ) = 5.33e11`
- 支配項は R1 の `[d u zeta*] Phi` で約 `7.58e75`

これは nonlinear formula の不整合ではなく、`rho` が zeta* の解析スケールより大きいことと、B.6 mixed structure を使わず standalone derivative bound を掛けたことを検出したものと解釈する。

## 成果物

CI は次を生成する。

- `results/reference/appendix_b_scale.csv`: B.13 から得る ordered `Lambda` scale
- `results/reference/appendix_b_brho.csv`: rho sweep、B.4/B.6 operator constants、multiplier norms、R1/R2 term-wise M、raw contraction proxy

## 未完了境界

次を満たすまでは Proposition B.2 の full contraction を完了扱いにしない。

1. `(1+T)^-1` の full Bρ operator bound を factorial degree-raising estimate から数値化する。
2. fixed-point map の post-`J2` / post-`J1` budget で B.6 mixed constants を各該当項に直接使用する。
3. `g=phi*/C` の complex-neighborhood bound を C の選択とともに閉じる。
4. radial/eta truncation を増やしたときの bound の安定性を確認する。
