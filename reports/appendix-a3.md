# Appendix A.3 moment closure 診断

このレポートは `include/appendix_a3_closure.hpp` の軽量数値回帰を GitHub Actions で再生成したもの。
論文の Appendix A.3 に対応する正規化された moment closure を検証する。完全な Theorem 4.6 profile assembly ではない。

## 主要結果

- `(A.19)` の `Kb = 0.245049619938`（論文中の評価 `0.20 < Kb <= 0.25` と整合）
- remainder を0とした principal amplitude: `1.01005026681`（論文の bracket `0.9 < Amp < 1.2` 内）
- `(A.15)` の `M,J` closure 最大残差: `1.735e-18`
- `(A.11)` の `I, pressure increment` closure 最大残差: `1.924e-63`
- `(A.16)` の `Qs(0)=Qp` 最大差: `2.045e-15`

## λ sweep

| lambda | Amp | c1(MJ) | c2(MJ) | c1(I/P) | c2(I/P) |
|---:|---:|---:|---:|---:|---:|
| 0.08 | 1.0100503 | -1.317e-02 | 5.605e-03 | -3.066e-33 | 3.120e-32 |
| 0.05 | 1.0100503 | -1.863e-06 | 7.498e-07 | -6.277e-39 | 5.665e-38 |
| 0.02 | 1.0100503 | 7.688e-21 | -3.018e-21 | -4.808e-50 | 3.848e-49 |

## 実装境界

現在は Appendix A.3 の有限次元 closure mechanism を式番号付きで再現している。
次段階では A.2 の全 radial schedule から `(A.14)` の pre-pulse discrepancies と `(A.19)` の remainder `E(Amp,eta)` を直接評価し、
この closure solver に入力する。これにより trial profile ではなく、outer profile 全体を通した `S(infinity)=0` の数値閉包へ進む。
