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

## 4. full schedule coupling と log-space 表現

A.2/A.3 の radial schedule から、second reserved patch と terminal heat tail の自然スケールを `include/appendix_a7_schedule_scales.hpp` で抽出する。ordered regime では両者が指数的に離れるため、`X_patch/X_tail` や `e_patch/e_tail` を通常の `double` として直接生成せず、すべて `log X`, `log e` で追跡する。

`lambda=1e-5`, `log h=-80` の基準では、second reserved patch と terminal tail の radial scale は `log X` で 10^6 オーダー離れる。このため、heat replacement により生じる実スケールの `Cp,S,I` discrepancy と、それを second reserved patch に換算した補正 target は通常の倍精度浮動小数点の下限を下回る。

この underflow は「補正が不要」という意味ではない。`include/appendix_a7_log_discrepancy.hpp` では (A.38) の leading expansion を使い、heat discrepancy を符号付き `log|Delta|` として保持する。cutoff 後の無限 tail は解析積分しており、小さい `h` でも有限窓打ち切り誤差を入れない。

さらに `include/appendix_a7_log_compensation.hpp` では、3つの target の最大 log magnitude を共通スケールとして取り出し、scaled Jacobian system を解く。係数そのものは符号付き log 表現で保持する。quadratic remainder は一次項よりさらに共通スケール `exp(L)` を1個余分に持つため、ordered regime では exact quadratic branch と linearized branch の差が表現限界よりはるかに小さいことを log bound で監査する。

`Cp` は axis pressure datum の increment そのものなので、同じ3-moment solve の第1成分を通じて pressure datum restoration も同じ scaled representation で追跡する。

## 5. 現在の実装境界

実装・CI回帰済み:

- (A.32) heat factor integral
- (A.34) derivative integral
- (A.35) derivatives at `Z=0`
- (A.37) heat-factor ODE regression
- (A.38) large-X expansion regression
- (A.39) terminal heat replacement ratio
- Proposition A.7 の second reserved patch 上の 3-bump exact quadratic moment map
- Lemma A.1 / A.2 を用いた `Cp,S,I` compensation solver
- terminal heat replacement の normalized discrepancy と second patch solver の runtime bridge
- A.2/A.3 full schedule から second patch / terminal tail の相対自然スケールを log-space で抽出
- ordered regime の heat discrepancy を signed-log で評価
- underflow-scale target を scaled Jacobian で解き、補正係数を signed-log で保持
- quadratic-to-linear remainder bound と pressure (`Cp`) target の整合性回帰

未完了:

- exact heat factor `H` による full discrepancy と large-X leading signed-log 式の誤差を、scaled asymptotic bound として end-to-end で閉じること
- replacement + compensation 後の `Cp,S,I` および axis pressure datum の recovery を、full schedule 全体の1つの結果型としてまとめること
- その end-to-end 診断を `results/reference/appendix_a6_a7_heat.csv` とレポートに固定し、Proposition A.7 を完了扱いにすること

したがって Proposition A.7 の full schedule coupling は大部分まで接続済みだが、まだ完了扱いにはしない。A.9 の残り3つの reserved intervals は引き続き未変更のまま保持する。
