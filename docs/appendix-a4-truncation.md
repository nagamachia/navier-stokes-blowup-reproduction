# Appendix A.4: radial truncation error audit

## 目的

A.4 pressure datum の誤差評価から、異なる刻み幅の解の差を使う Richardson estimator を外す。

現在の実装は `appendix_a4_pressure_enclosure` で pressure datum を **1回だけ**評価し、次の a-priori bound を使う。

- constant-`l` 区間: 解析積分なので truncation error なし
- transition 区間: RK4 の global `O(h^4)` error を 5階微分 majorant から評価
- A.10: composite Simpson の `O(h^4)` error を integrand の4階微分 majorantから評価

## transition の微分 bound

`g=E^2/P_*^2` とし、`r=(log g)'=2l-1` とすると

`g' = r g`。

したがって

`g^(5)/g = r^5 + 10 r^3 r' + 15 r (r')^2 + 10 r^2 r'' + 10 r' r'' + 5 r r''' + r''''`。

A.4 の transition はすべて smooth step `s` の affine combination なので、`s^(k)` の上界から `r,...,r''''` の上界へ直接移せる。

現在の reproduction majorant は

- `|s'| <= 9`
- `|s''| <= 100`
- `|s'''| <= 4000`
- `|s''''| <= 120000`
- `|s'''''| <= 7000000`

を大きめに丸めた real-axis lemma として分離している。これらは**原論文の定数ではない**。

この lemma を使うと transition 全長7に対して、sharp RK4 principal coefficient を使わずさらに保守的に

`E_RK4 <= 7 h^4 G M5`

と置く。ここで `G=1.25`、`M5` は上の Bell polynomial majorant である。

## A.10 Simpson

A.10 の integrand は

`g(y,eta)=c(y) exp(2(theta(y)-1)a(eta))`

の形で、`theta(y)=1-s(y/Tf)`。

`y` 微分のたびに smooth-step derivative は `Tf^(-k)` を伴う。現在は transition の `M5` を一段強い共通 majorant として流用し、標準 composite Simpson bound

`E_S <= Tf h^4 max|F''''| / 180`

を適用する。

## eta 微分への伝播

A.10 前の eta dependence は `f(eta)^2` だけなので、radial integral error は明示的な `J_eta`, `J_etaeta` の係数を通して `Pi0_eta`, `Pi0_etaeta` に伝播する。

`|eta|<=1` で

- `|J_eta|<=1`
- `|J_etaeta|<=2`
- `|-2J_eta|<=2`
- `|4J_eta^2-2J_etaeta|<=8`

を使っている。A.10 の eta multiplier にも同じ保守 bound を使う。

## Richardson から変わった点

旧実装は `H`, `H/2`, `H/4` の3解を計算し、差を `/15` して safety factor を掛けていた。

新実装は center accuracy を維持するため旧 `H/4` と同じ刻みを使うが、`H` と `H/2` は計算せず、resolution difference は error estimate に一切使わない。

## 現在の rigor 境界

truncation-error **式**は a-priori derivative majorant に置き換わった。一方、上記5個の smooth-step derivative 定数は現在 reproduction lemma としてコードに固定した段階であり、まだ interval jet / symbolic interval arithmetic による独立証明を付けていない。

したがって現時点では「Richardson 依存を除き、解析的 truncation bound に必要な微分上界へ問題を還元した」が正確な位置づけである。次の rigor step は smooth step の1〜5階微分を区間演算で全 `t in [0,1]` certification し、この lemma の仮定依存を外すことである。
