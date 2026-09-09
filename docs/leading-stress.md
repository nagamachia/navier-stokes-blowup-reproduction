# Leading stress の論文抽出メモ

この文書は OpenAI 論文 Section 4.1 の式 (4.8)–(4.11) を、後続の数値実装用に日本語で整理したもの。ここに書く式は論文PDFを直接確認したもの。

前提として式 (4.7) の `E,U,Pi,V0` が成立し、`phi,U,Pi` は軸で滑らか、`phi > 0` とする。

## 式 (4.8): 補助量

論文は

`H = sqrt(2X) E`

`F = E / sqrt(2X)`

`l = D_X log H`

を置く。

さらに

`W = 1 - 2 D eta A_X(U) - d partial_eta A_X(U)`

`H_c = D eta + d U`

を定義する。

ここで

`D_X = X partial_X`

`d = 1 - eta^2`

`D = 1/2 - h`

である。

`H` は angular momentum `r u_theta^(0)` の profile。`F` は `E` から軸上の `sqrt(2X)` 因子を除いた量で、`l` は `H` の対数半径勾配。

## 式 (4.9): stress source

2つの scalar source `S_q`, `S_n` を

`D_X Q_s + (1+l) Q_s = S_q`

`D_X N_s + N_s = S_n`

の右辺として定義する。

`S_q = -W l - h(1 - 2 eta U) - H_c partial_eta(log E)`

`S_n = -W D_X U - A(1 - 2 eta U) U - H_c partial_eta U`

`      - d partial_eta Pi + 4 A eta Pi + 2 eta D_X Pi`

ここで `A = 1/2 + h`。

## 式 (4.10): radial integration

一般解は積分定数 `C_Q(eta)`, `C_N(eta)` を含み

`Q_s = [C_Q(eta) + integral_0^X H(x,eta) S_q(x,eta) dx] / [X H(X,eta)]`

`N_s = [C_N(eta) + integral_0^X S_n(x,eta) dx] / X`

となる。

軸で滑らかに延長するため、論文は

`C_Q = C_N = 0`

を選ぶ。変数 `x=sX` を使うと

`Q_s = F(X,eta)^(-1) integral_0^1 s F(sX,eta) S_q(sX,eta) ds`

`N_s = integral_0^1 S_n(sX,eta) ds`

となる。

軸上極限は

`Q_s(0,eta) = S_q(0,eta)/2`

`N_s(0,eta) = S_n(0,eta)`

である。

## 式 (4.11): 2成分 stress coefficient

論文は `(r theta, r z)` 成分の2成分係数を

`T_0 = F (p_s - s_bar)`

と置く。

`p_s = (X Q_s / L, X N_s / (L E))`

`s_bar = (a, -b_s)`

`a = 1 - 2 D_X log E = 2 - 2 l`

`b_s = 2 D_X U / E`

である。

ここで `p_s` は積分された inviscid contribution、`-F s_bar` は radial-viscosity contribution。物理 stress は

`T = q^(-A-1/2) T_0`

として現れる。

## 数値実装時の注意

1. `E ~ sqrt(2X) F` なので、`E` や `log E` を直接 `X=0` で評価してはいけない。
2. 論文自身が `(log E)_eta = F_eta/F` を使い、軸で滑らかに延長している。
3. `Q_s`, `N_s` は式 (4.10) の `s in [0,1]` 表現を使うと `X=0` を含めて安定に実装しやすい。
4. exact `E,U,Pi` は Theorem 4.6 と Appendices A/B の構成から得る必要がある。テスト用 manufactured profile と exact construction を区別する。
5. この `T_0` が Sections 6–7 の oscillatory waves で実現すべき covariance stress の対象になる。

## 次の実装順

- generic profile に対する `H,F,l,W,H_c`
- `S_q,S_n`
- 式 (4.10) の smooth radial integral
- `Q_s,N_s,p_s,a,b_s,T_0`
- manufactured profile で ODE identity と軸上極限を回帰テスト
- その後に Theorem 4.6 の actual profile construction を数値化する
