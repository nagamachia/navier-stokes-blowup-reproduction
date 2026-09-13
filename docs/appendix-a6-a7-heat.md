# Appendix A.6–A.7 heat exterior 数値実装ノート

原典: *FINITE TIME BLOWUP FOR NAVIER–STOKES* (OpenAI, 2026), Appendix A.6–A.7, Lemma A.6 / Proposition A.7。

## 1. Lemma A.6: exact heat factor

`A = 1/2 + h`, `d = 1 - eta^2` とし、原典 (A.32) の heat factor

`H(Z) = Gamma(1+h)^(-1) integral_0^infinity exp(-v) v^h (1+Zv)^(-h) dv`

を `include/appendix_a6_heat.hpp` に実装した。

弱い endpoint power `v^h` をそのまま等間隔 Simpson 積分すると小さい `h` で収束が遅いため、数値積分では `v=u^8` と変数変換し、原点を平滑化してから composite Simpson quadrature を適用する。

同じ積分表示から (A.34) の `m` 階微分も直接評価し、(A.35)

`H^(m)(0) = (-1)^m (h)_m (1+h)_m`

と照合する。

さらに (A.37)

`Z^2 H'' + (1 + 2(1+h)Z) H' + h(1+h) H = 0`

の residual、`H>0`, `H'<0`, `0 <= -ZH'/H < h` を回帰する。

## 2. heat exterior ratio と terminal replacement

power tail に対する exact heat profile の比は

`E_heat / E_pow = H(2d/X)`

である。`appendix_a6_heat_ratio()` がこれを評価する。

(A.38) の large-X expansion

`H(2d/X) = 1 - 2h(1+h)d/X + O(X^-2)`

は、一次近似との差に `X^2` を掛けた量が bounded / asymptotically stable であることをテストする。

(A.39) の terminal replacement は local coordinate `y=log(X/X_tail)` で、原典と同じく `y<=0.2` では変更せず、`y>=0.5` で full heat ratio に到達する smooth cutoff として実装した。

## 3. Proposition A.7: second reserved patch の 3-moment compensation

heat replacement は `U=0` のままなので、式 (4.15) の `M` と `J` は構造的に変化しない。一方で、次の3量を戻す必要がある。

- pressure increment `Cp = integral E^2/(2X) dX`
- `S = integral (U^2 - E^2/2) dX`。ここでは `U=0`
- angular moment `I = integral sqrt(2X) E dX`

A.9 の **second reserved patch**（heat compensation 用、log-width 5）を normalized coordinate `x=X/X_*` で `[1,e^5]` とし、互いに分離した3つの nonnegative smooth bumps `beta_j` を置く。

基準 profile は intermediate power law

`E0/e_* = f(eta) x^(-1/2-lambda)`, `f(eta)=1/(1+eta^2)`

で、補正を

`delta E/e_* = sum_j c_j beta_j(x)`

とする。

`AppendixA7HeatCompensation` は `c -> (Delta Cp, Delta S, Delta I)` の **exact degree-two map** を積分し、その zero-state Jacobian と quadratic remainder を構成する。線形化した3行の重みはそれぞれ

- `f x^(-3/2-lambda)`
- `-f x^(-1/2-lambda)`
- `sqrt(2x)`

となり、相異なる power weight と分離 support により Lemma A.1 の可逆性機構を数値確認できる。

係数は既存の Lemma A.2 実装 `solve_small_quadratic_moment_system()` を用いて zero-start branch から解く。回帰では `eta` sweep に対して normalized Jacobian の非退化、既知の小係数から生成した3 discrepancy の exact recovery、`f(eta)` に由来する `eta -> -eta` 対称性を確認する。

## 4. 現在の実装境界

実装済み:

- (A.32) heat factor integral
- (A.34) derivative integral
- (A.35) derivatives at `Z=0`
- (A.37) heat-factor ODE regression
- (A.38) large-X expansion regression
- (A.39) terminal heat replacement ratio
- Proposition A.7 の second reserved patch 上の 3-bump exact quadratic moment map
- Lemma A.1 / A.2 を用いた `Cp,S,I` compensation solver

未接続:

- A.2/A.3 の full radial schedule から terminal tail の実スケール `X_K`, `e_K` を取り出す処理
- heat replacement が tail 全体に与える `(Delta Cp, Delta S, Delta I)` を full schedule の単位系で積分する処理
- その実 discrepancy を second reserved patch solver へ渡して、replacement + compensation 後の global moments と axis pressure datum が元値へ戻ることの end-to-end 回帰

したがって現時点では Proposition A.7 の有限次元 compensation mechanism までを実装済みとし、full schedule を含む Proposition A.7 全体の完了とは扱わない。
