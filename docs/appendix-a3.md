# Appendix A.3 数値実装ノート

原典: `FINITE TIME BLOWUP FOR NAVIER–STOKES` (OpenAI, 2026), Appendix A.3 `Closing the moment integrals`。

この実装では、Appendix A.3 の有限次元 closure mechanism を、全 radial schedule から切り離して数値回帰できる形にしている。現時点では Theorem 4.6 の完全 profile assembly ではない。

## 1. axial pulse と M,J closure

A.2 の pulse は

`xi_b = lambda y`

`phi_b(xi) = integral_0^xi sigma(v/0.02) dv`

`R0(xi) = phi_b(xi) [1-sigma(xi-10)]`

`Rb = Amp(eta) R0 + c1 beta1 + c2 beta2`

`U = E Rb`

である。`beta1,beta2` は幅 0.3、中心 `13/lambda-3`, `13/lambda-1` の同一translate。

Appendix A.3 の (A.15) では、M,J をそれぞれ first bump center で正規化すると重みは

`s1 = 1/2-lambda`

`s2 = 1/2-2lambda`

`exp(s1(y-y1)), exp(s2(y-y1))`

となる。実装ではこの正規化を直接用い、巨大な `exp(O(1/lambda))` を生成しない。

`appendix_a3_close_pulse_MJ()` は trial amplitude と pre-pulse discrepancy から `c1,c2` を2×2線形系で求め、補正後の `M=J=0` を回帰する。

## 2. angular I / pressure correction

(A.11) の relative E correction は、uniform `l=-lambda` 区間で2つのE bumpを使う。

正規化されたI rowの重みは `exp((1-lambda)(y-y1))`、pressure rowの線形化は `exp((-1-2lambda)(y-y1))`。

pressure increment はEの二乗に依存するため、実際の方程式は二次である。`appendix_a3_close_angular_I_pressure()` は Lemma A.2 の zero-start fixed-point branch をそのまま使い、Iを指定しつつ total pressure incrementを保持する。

## 3. amplitude equation

(A.19) の主要方程式は

`Amp^2 Kb - (1-exp(-26))/4 + Error(Amp,eta) = 0`

`Kb = integral_0^13 exp(-2 xi) R0(xi)^2 dxi`

である。原典は `0.20 < Kb <= 0.25` を示し、十分小さいlambdaで `Amp in (0.9,1.2)` に一意なrootを取る。

`appendix_a3_Kb()` と `appendix_a3_solve_amplitude()` はこの式を直接数値化する。現段階では full radial schedule から `Error(Amp,eta)` を再構成せず、callbackとして与える。remainder=0 の principal equation と小さなeta依存remainderを回帰する。

## 4. terminal Qs

(A.16) のterminal representation

`Qs(y) = integral_y^3 exp(integral_y^v (1+l(s)) ds) [fo'(v)/fo(v)] dv`

を `appendix_a3_terminal_Qs()` で評価する。(A.13) の `Qp` と `Qs(0)=Qp`、endpointの `Qs(3)=0` を数値照合する。

また (A.17) の `l=-1` hold interval `4 log(1/h)` に対する exact factors

`(XE^2)_end/(XE^2)_start = h^8`

`E_end/E_start = h^6`

も回帰する。

## 5. 現在の境界と次の実装

ここまでで A.3 の「有限次元 closure の機構」は実装済み。一方、faithful outer-profile reproduction には次がまだ必要。

- A.2 の全stageを単一のglobal log-radial scheduleとして接続
- そのscheduleから (A.14) の pre-pulse discrepanciesを直接積分
- profile interpolation / exterior releaseまで含めて (A.19) の `Error(Amp,eta)` を直接評価
- 得られた `Amp(eta)`、pulse correction、angular correctionをprofileへ実際に反映
- 全radial integralで (A.8) の2条件を再確認

次はこの統合を行い、その後 A.4 の pressure datum `Pi0(eta)` へ進む。
