# Appendix A.4 pressure datum 数値実装ノート

原典: `FINITE TIME BLOWUP FOR NAVIER–STOKES` (OpenAI, 2026), Appendix A.4 / Lemma A.5, 式 (A.21)–(A.23)。

## 1. pressure datum

`y = log(X/X_R)` とすると、(A.21) の axis pressure datum は

`Pi0(eta) = -1/2 integral_{-infinity}^{infinity} E(y,eta)^2 dy`

である。`dX/X = dy` なので `X_R` は積分から消え、`Pi0` は radial scale の選び方に依存しない。

実装 `appendix_a4_pressure_datum()` は `Pi0/P*^2` を基本量にして積分する。これにより proof hierarchy で大きく取る `P*` を最後まで直接二乗せず、必要な場合だけ `appendix_a4_scale_pressure()` で物理スケールへ戻す。

## 2. A.2/A.3 schedule との接続

pressure 積分には既存の A.2 schedule をそのまま使う。

- `y<=0`: (A.7) reference inner branch
- first slope transition
- axial reduction interval `Td`
- intermediate `l=-lambda` への transition と `Tw`
- axial pulse interval `13/lambda`
- (A.10) eta-dependence removal
- (A.11) angular adjustment interval
- exterior `-lambda -> -1`, `l=-1` hold, `-1 -> -h`
- `l=-h` hold
- (A.12) terminal collar
- terminal power-law tail to `y=infinity`

A.11 の2つの relative-E bumps は total pressure increment を保存するよう選ばれるため、(A.21) の datum 評価では bumps 自体を省略できる。ただし A.11 の radial interval は削除しない。

axial pulse の `U`, `Amp(eta)`, M/J correction は `E` を変更しないため pressure datum の integrand には直接入らない。

## 3. (A.22) の意味

inner branch は

`E/P* = f(eta) exp(y/10)`, `f(eta)=1/(1+eta^2)`

なので

`integral_{-infinity}^0 E^2/P*^2 dy = 5 f(eta)^2`。

従って inner branch だけで

`Pi_inner/P*^2 = -(5/2) f(eta)^2`

となる。外側の integrand も非負なので full datum は

`Pi0(eta) <= -(5/2) P*^2 f(eta)^2`

を満たす。これは (A.22) の bound であり、full `Pi0` を右辺と同一視してはいけない。

同時に full schedule は eta について偶関数で、原典の符号条件

`eta Pi0'(eta) > 0  (eta != 0)`

を数値回帰する。

## 4. eta derivatives

(A.10) までは `E^2 = c(y)^2 f(eta)^2`。A.10 中は

`E = c(y) f(eta)^theta(y)`

と書ける。`J=log(1+eta^2)` とすると

`partial_eta E^2 = -2 theta J_eta E^2`

`partial_eta^2 E^2 = (4 theta^2 J_eta^2 - 2 theta J_etaeta) E^2`

を integrand と同時に積分する。A.10 終了後は eta dependence が消えるため、後続stageは `Pi0_eta`, `Pi0_etaeta` に寄与しない。

テストでは `Pi0_eta`, `Pi0_etaeta` を finite difference と照合する。この derivative は Appendix B の `Z*` に直接必要になる。

## 5. (A.23) inner pressure

(A.23) は finite radius で

`Pi(X,eta) = -1/2 integral_y^infinity E(v,eta)^2 dv`, `y=log(X/X_R)`

を与える。reference inner branch `y<=0` では厳密に

`Pi(y,eta)/P*^2 = Pi0(eta)/P*^2 + (5/2) f(eta)^2 exp(y/5)`

となるので `appendix_a4_inner_pressure_normalized()` として実装した。

`y -> -infinity` で `Pi -> Pi0`、また `partial_y Pi = E^2/2 > 0` である。

## 6. 回帰と境界

`tests/appendix_a4_pressure.cpp` は以下を検証する。

- `Pi0` の偶対称性
- `Pi0_eta` の奇対称性
- (A.22) の inner bound
- `eta Pi0_eta > 0`
- analytic eta derivatives と finite difference の一致
- (A.23) inner pressure の極限・単調性
- `P*^2` scaling の有限性

CIでは `results/reference/appendix_a4_pressure.csv` と `reports/appendix-a4-pressure.md` を生成する。

これは有限パラメータでの式・schedule回帰であり、Lemma A.5 の存在証明や theorem hierarchy 全体を数値的に証明するものではない。

次は `Pi0`, `Pi0_eta` を Appendix B.1 の `Z*` に接続し、analytic axis profile と moment matching を実装する。
