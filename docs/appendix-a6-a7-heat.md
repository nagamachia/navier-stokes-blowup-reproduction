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

moderate な normalized regime では、既存の Lemma A.2 実装 `solve_small_quadratic_moment_system()` を用いて3 discrepancy の exact numerical recovery を回帰している。

## 4. full schedule coupling と log-space 階層

A.2/A.3 の radial schedule から、second reserved patch と terminal heat tail の自然スケールを `include/appendix_a7_schedule_scales.hpp` で抽出する。ordered regime では両者が指数的に離れるため、`X_patch/X_tail` や `e_patch/e_tail` を通常の `double` として直接生成せず、すべて `log X`, `log e` で追跡する。

`lambda=1e-5`, `log h=-80` の基準では、second reserved patch と terminal tail の radial scale は `log X` で 10^6 オーダー離れる。このため、heat replacement により生じる実スケールの `Cp,S,I` discrepancy と、それを second reserved patch に換算した補正 target は通常の倍精度浮動小数点の下限を下回る。

この underflow は「補正が不要」という意味ではない。`include/appendix_a7_log_discrepancy.hpp` では (A.38) の leading expansion を使い、heat discrepancy を符号付き `log|Delta|` として保持する。cutoff 後の無限 tail は解析積分しており、小さい `h` でも有限窓打ち切り誤差を入れない。

さらに重要なのは、ordered regime では `Cp,S,I` の patch target 同士にも 10^6 オーダーの log-scale 差が生じることである。最大の target を O(1) に正規化すると、最小の target は倍精度では消える。また、最大 target に対する quadratic term は十分小さい一方で、その quadratic term 自体が最小 target より大きくなり得る。このため、通常の倍精度3×3 solve を「full ordered-regime exact compensation」と呼ぶことはできない。

`include/appendix_a7_log_compensation.hpp` はこの事実を隠さず、**scale audit** として実装する。最大 target に対する dominant linear direction、target の log spread、quadratic-to-dominant ratio、quadratic absolute log magnitude を計算し、全3成分を倍精度で同時解像可能かを明示する。ordered 基準では `all_components_resolvable_in_double=false` になる。

したがって現時点の数値的主張は次の2段階に分ける。

1. moderate normalized regime: exact degree-two map と3-moment compensation を end-to-end で数値回帰できる。
2. paper の ordered regime: full schedule の指数的階層と dominant correction / nonlinear scale を log-space で監査できるが、全3 target の componentwise exact recovery を倍精度だけで直接計算したとは主張しない。

`Cp` は axis pressure datum の increment そのものなので、pressure restoration も同じ精度階層の制約を受ける。

## 5. 現在の実装境界

実装・回帰済み:

- (A.32) heat factor integral
- (A.34) derivative integral
- (A.35) derivatives at `Z=0`
- (A.37) heat-factor ODE regression
- (A.38) large-X expansion regression
- (A.39) terminal heat replacement ratio
- Proposition A.7 の second reserved patch 上の 3-bump exact quadratic moment map
- moderate normalized regime の `Cp,S,I` exact compensation solver
- terminal heat replacement の normalized discrepancy と second patch solver の runtime bridge
- A.2/A.3 full schedule から second patch / terminal tail の相対自然スケールを log-space で抽出
- ordered regime の heat discrepancy を signed-log で評価
- ordered regime の target log hierarchy、dominant linear correction、quadratic scale の監査
- 診断CSVへの ordered target logs / scale spread / nonlinear bound の出力

未完了:

- exact heat factor `H` による full discrepancy と large-X leading signed-log 式の誤差を scaled asymptotic bound として閉じること
- ordered regime の階層補正について、解析的な contraction / hierarchy argument と数値監査を接続し、`Cp,S,I` と axis pressure datum restoration をどこまで「再現済み」と判定できるかを明文化すること
- その上で Proposition A.7 を完了扱いにすること

したがって Proposition A.7 の full schedule coupling は進行中であり、まだ完了扱いにはしない。A.9 の残り3つの reserved intervals は引き続き未変更のまま保持する。
