# 実装進捗

この文書は、論文のどこまでを「式番号付きで抽出したか」「コード化したか」「数値検証したか」を区別して記録する。

## 完了

### 数値ソルバ基盤

- Fourier 微分
- 非圧縮射影
- 2/3 dealiasing
- Fourier モードの粘性減衰
- 直接 DFT の Taylor–Green 回帰
- FFTW 版 Taylor–Green 実行系
- CSV 診断量
- GitHub Actions による自動回帰

### 論文構成

- 式 (3.2)/(4.1): 類似座標 `q, eta, X`
- Lemma 4.1, 式 (4.2): `T_b`, `Z_b`
- 式 (4.6)–(4.7): `A_X`, `V0`, `Pi`
- 式 (4.8)–(4.11): leading stress
- 式 (4.4)–(4.5): Cartesian 軸正則性
- 式 (4.15): five cumulative radial moments `M,I,J,S,Cp`
- Appendix A.1, Lemma A.1: finite moment matrix / invertibility
- Appendix A.1, Lemma A.2: quadratic moment correction
- Corollary A.3: 2つの `U` bumps + 3つの `E` bumps による five-moment exact quadratic map
- Appendix A.2 staged outer profile の radial primitives
  - (A.5) smooth step `sigma`
  - (A.6) parameter hierarchy を `log P*`, `log h` で安全に表現
  - (A.7) reference inner power law
  - axial reduction / (A.9) intermediate power law の基準stage
  - (A.9) に置く4つの reserved correction interval の位置
  - (A.10) profile interpolation
  - (A.12) terminal factor / slope
  - (A.13) `Qp`
- Appendix A.3 finite-dimensional moment closure mechanism
  - axial pulse `phi_b`, `R0`
  - (A.15) の2-bump `M=J=0` closure
  - (A.11) の2-bump `I` / pressure increment closure
  - (A.19) の `Kb` と amplitude root solver
  - (A.16) terminal `Qs` representation
  - (A.17) の `h^8`, `h^6` exact factors
  - ordered asymptotic regime `lambda=1e-5` で full schedule の A.19 bracket/root を回収
  - `Amp(eta)` を eta sweep で解き、global closure を CI 回帰
- Appendix A.2 → A.3 の基準schedule積分
  - temporary reference power law → first transition → axial reduction → (A.9) hold を正規化ODEで接続
  - `m=M/(XE)`, `j=J/(XHE)`, `s=S/(XE^2)`, `rI=I/(XH)` を直接発展させる実装
  - raw `X,M,J,S` を生成せず巨大な対数半径を扱える
  - scheduleから得た `M,J` discrepancy を (A.15) two-bump closure に接続
- Appendix A.4 / Lemma A.5 pressure datum
  - (A.21) `Pi0(eta) = -1/2 integral E(y,eta)^2 dy` を A.2/A.3 の全 radial schedule から評価
  - `Pi0/P*^2` を基本量とする log-space 積分で巨大な `P*` を安全に扱う
  - (A.10) 中の `Pi0_eta`, `Pi0_etaeta` を integrand の解析微分から同時積分
  - (A.22) の偶対称性、`eta Pi0'(eta)>0`、inner-branch bound を回帰
  - (A.23) の `y<=0` inner pressure を厳密形で実装
  - A.11 の pressure-preserving angular bumps は datum 積分では省略するが、A.11 区間そのものは保持
- Appendix A.6 / Lemma A.6 exact heat exterior primitives
  - (A.32) `H(Z)` の積分表示を実装
  - `v=u^8` 変数変換で `v^h` の endpoint quadrature を安定化
  - (A.34) の `H^(m)`、(A.35) の `Z=0` derivatives を数値照合
  - (A.37) ODE residual、`H>0`, `H'<0`, `0<=-ZH'/H<h` を回帰
  - (A.38) `E_heat/E_pow = H(2d/X)` の `O(X^-2)` remainder を回帰
  - (A.39) の terminal heat replacement ratio を実装
- Proposition A.7 heat compensation の有限次元コア
  - A.9 の second reserved patch のみを normalized `x in [1,e^5]` として使用
  - 3つの separated nonnegative `E` bumps に対する exact quadratic `(Delta Cp,Delta S,Delta I)` map
  - zero-state Jacobian の3つの異なる power weights と非退化を eta sweep で確認
  - Lemma A.2 zero-start solver で synthetic small discrepancy を machine precision 近くまで回収
  - `U=0` のため `M,J` は構造的に不変
  - `eta -> -eta` 対称性を回帰
- Proposition A.7 full schedule coupling の数値部品
  - terminal heat replacement の normalized discrepancy → second patch solver bridge を runtime 回帰
  - bridge 内の `appendix_a6_terminal_replacement_ratio` 引数順誤りを発見・修正
  - A.2/A.3 full schedule から second reserved patch / terminal tail の `log X`, `log e` 自然スケールを抽出
  - ordered regime (`lambda=1e-5`, `log h=-80`) で patch と tail が `log X` で 10^6 オーダー離れることを確認
  - heat discrepancy と patch target を signed-log 表現で保持
  - target 3成分の log spread、dominant linear correction、quadratic-to-dominant ratio を監査
  - `all_components_resolvable_in_double=false` を明示し、倍精度での componentwise exact recovery を誤って主張しない
  - heat factor Taylor remainder から leading discrepancy の相対誤差 bound を導出し、ordered regime では `log(relative error)<-10^6` を回帰

### 数式監査・数値表現で検出した修正

- 類似座標の消去形は `q - z^2 q^(2h) = tau` と確認し、初期転記を修正した。
- Corollary A.3 のテストで power weights が重複する退化 exponent を避けるよう修正した。
- Appendix A.2 の巨大パラメータは `log P*`, `log h` を基本表現にした。
- Appendix A.3 の pulse moments は first bump center で正規化し、`exp(O(1/lambda))` を直接生成しない。
- A.11 の relative-E correction は局所 moment jump と global slope evolution を二重計上しないよう分離した。
- **重要な訂正:** (A.9) 後に列挙される4区間は A.3 closure のために使う補正ではない。Proposition A.4 はこれらを unaltered のまま残して (A.8) を成立させ、その後の cone realization / heat compensation / higher-order background / phase-averaged mean corrections が使用する。
- `lambda=0.05` を用いた初期 global test では `S/(XE^2)` が巨大化して A.19 bracket が消えた。これは reserved patch 不足ではなく、`Tw=60 log(1/lambda)` を含む ordered asymptotic construction に対して lambda が十分小さくなかったため。A.14 で使う `exp(-2 lambda Tw)=exp(-120 lambda log(1/lambda))` は lambda -> 0 で 1 に近づくが、lambda=0.05 では約 1.6e-8 まで落ちる。
- (A.22) の `-(5/2)P*^2 f(eta)^2` は full pressure datum の等式ではなく、`y<=0` reference inner branch から来る上界として扱う。
- (A.32) の `v^h` endpoint は naive uniform-v Simpson では小さい `h` の収束が遅いため、`v=u^8` で平滑化してから積分する。
- A.7 full schedule の ordered regime では `Cp,S,I` patch target のスケール差が 10^6 オーダーに達する。最大 target を O(1) に正規化すると最小 target は double で消え、dominant correction の quadratic term が最小 target より大きくなり得る。したがって ordered regime の full exact compensation は単純な倍精度3×3 solve では検証できず、階層的 nonlinear correction / contraction argument の接続が必要。

### 成果物パイプライン

main 更新時に軽量基準計算を実行し、Taylor–Green検証、3D時系列VTK、可視化、レポートを GitHub に残す。

- Appendix A.3: `results/reference/appendix_a3_closure.csv`, `results/reference/appendix_a3_global.csv`, `reports/appendix-a3.md`
- Appendix A.4: `results/reference/appendix_a4_pressure.csv`, `reports/appendix-a4-pressure.md`
- Appendix A.6/A.7: `results/reference/appendix_a6_a7_heat.csv`, `docs/appendix-a6-a7-heat.md`

A.6/A.7 CSV は `eta` sweep に対して heat ratio、(A.37) ODE residual、normalized 3x3 Jacobian determinant、moderate-regime moment recovery error に加え、ordered regime の target signed-log、target log spread、quadratic scale、leading-asymptotic relative-error bound、double で全成分を解像可能かを保存する。

## 現在地: Proposition A.7 full schedule coupling

Lemma A.6 の exact heat factor / terminal replacement、second reserved patch 上の3-moment finite-dimensional compensation mechanism、full schedule scale extraction、ordered-regime signed-log discrepancy と leading-asymptotic error bound まで実装・回帰済み。

moderate normalized regime では replacement discrepancy → second patch compensation の end-to-end recovery を通常の倍精度で回帰できる。一方、論文の ordered regime では3 target のスケール差が極端で、full componentwise exact solve は倍精度の表現能力を超えることが数値的に確認された。

したがって Proposition A.7 全体はまだ完了扱いにしない。次に必要なのは、ordered regime の dominant correction が作る quadratic remainder をさらに小さい補正階層へ送る **hierarchical nonlinear correction / contraction** を、Lemma A.2 の理論と対応する監査量として実装すること。その上で `Cp,S,I` と axis pressure datum restoration のうち、数値で直接確認する部分と解析 bound で保証する部分を1つの end-to-end 結果型にまとめる。

A.9 の残り3つの reserved intervals は引き続き未変更のまま保持する。

次の順序:

1. **Appendix A.1 moment correction primitives — 実装済み**
2. **式 (4.15) cumulative moments — 実装済み**
3. **Corollary A.3 five-bump exact quadratic moment map — 実装済み**
4. **Appendix A.2 radial primitives / (A.10)–(A.13) — 実装済み**
5. **Appendix A.3 finite-dimensional closure mechanism — 実装済み**
6. **A.2 基準scheduleの正規化ODE — 実装済み**
7. **A.9 reserved correction intervals — 位置実装済み、second patch のみ A.7 compensation core で使用開始**
8. **ordered asymptotic regime で corrected A.14 / global A.19 bracket — 数値確認済み**
9. **A.3 `Amp(eta)` global root / eta sweep — 実装・CI回帰済み**
10. **Appendix A.4 `Pi0(eta)` / (A.21)–(A.23) — 実装・CI回帰済み**
11. **Lemma A.6 exact heat factor / (A.32)–(A.39) — 実装・CI回帰済み**
12. **Proposition A.7 3-moment finite-dimensional compensation core — 実装・CI回帰済み**
13. **Proposition A.7 full schedule scale / signed-log / asymptotic-error audit — 実装・CI回帰済み**
14. Proposition A.7 hierarchical nonlinear correction / contraction audit
15. Appendix B analytic axis profile と moment matching
16. Appendix C admissible stress cone realization
17. Appendix C.2 後に first reserved patch で five moments を復元
18. higher-order / mean corrections で third/fourth reserved patches を使用
19. Theorem 4.6 profile assembly

任意の surrogate profile を論文構成そのものとして扱わない。

## その後

1. Sections 6–7 の oscillatory pulse
2. Propositions 7.5–7.6 の stress realization
3. 残差改善
4. Section 10 の localization と smooth forcing
5. 必要に応じて 64^3 → 128^3 → 256^3 解像度研究
