# Appendix A.3 moment closure 診断

このレポートは Appendix A.2 の staged radial profile を A.3 の closure solver に接続した軽量数値回帰を GitHub Actions で再生成したもの。
`M/(XE)`, `J/(XHE)`, `S/(XE^2)` を直接発展させるため、巨大な `X` や `P*` を明示的に生成せず theorem-scale の対数半径を扱える。
完全な Theorem 4.6 profile assembly ではなく、現在は axial pulse 開始までの global schedule と A.15 closure を直接接続している。

## 主要結果

- `(A.19)` の `Kb = 0.245049619938`（論文中の評価 `0.20 < Kb <= 0.25` と整合）
- remainder を0とした principal amplitude: `1.01005026681`（論文の bracket `0.9 < Amp < 1.2` 内）
- A.2 schedule から A.15 の第1 bump 中心へ運んだ pre-M 最大値: `4.870e-71`
- A.2 schedule から A.15 の第1 bump 中心へ運んだ pre-J 最大値: `3.123e-60`
- `(A.15)` の `M,J` closure 最大残差: `1.735e-18`
- `(A.11)` の `I, pressure increment` closure 最大残差: `1.924e-63`
- `(A.16)` の `Qs(0)=Qp` 最大差: `2.045e-15`

## λ・η sweep

| lambda | eta | m(Xp) | j(Xp) | s(Xp) | pre-M | pre-J | c1(MJ) | c2(MJ) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.08 | -0.5 | -6.039e-42 | -1.113e-36 | -4.678e+11 | -4.870e-71 | -3.123e-60 | -1.321e-02 | 5.623e-03 |
| 0.08 | 0 | 0.000e+00 | 0.000e+00 | -4.678e+11 | 0.000e+00 | 0.000e+00 | -1.321e-02 | 5.623e-03 |
| 0.08 | 0.5 | 6.039e-42 | 1.113e-36 | -4.678e+11 | 4.870e-71 | 3.123e-60 | -1.321e-02 | 5.623e-03 |
| 0.05 | -0.5 | -1.944e-49 | -1.533e-45 | -9.830e+08 | -1.155e-99 | -3.468e-90 | -1.898e-06 | 7.648e-07 |
| 0.05 | 0 | 0.000e+00 | 0.000e+00 | -9.830e+08 | 0.000e+00 | 0.000e+00 | -1.898e-06 | 7.648e-07 |
| 0.05 | 0.5 | 1.944e-49 | 1.533e-45 | -9.830e+08 | 1.155e-99 | 3.468e-90 | -1.898e-06 | 7.648e-07 |
| 0.02 | -0.5 | -3.018e-63 | -3.204e-61 | -2.696e+05 | -4.029e-198 | -1.782e-190 | -4.808e-22 | 1.834e-22 |
| 0.02 | 0 | 0.000e+00 | 0.000e+00 | -2.696e+05 | 0.000e+00 | 0.000e+00 | -4.808e-22 | 1.834e-22 |
| 0.02 | 0.5 | 3.018e-63 | 3.204e-61 | -2.696e+05 | 4.029e-198 | 1.782e-190 | -4.808e-22 | 1.834e-22 |

## 実装境界

A.2 の reference power law → first transition → axial reduction → intermediate entry → `(A.9)` hold を一つの正規化ODEとして積分し、
得られた `(A.14)` の pre-pulse discrepancies を `(A.15)` の two-bump solver へ直接入力するところまで統合済み。
次は pulse 本体・A.10 profile interpolation・A.11 angular correction・exterior transition を同じ正規化状態に追加し、
`S(infinity)` を trial amplitude ごとに直接評価して `(A.19)` の remainder を callback ではなく global schedule から構成する。
