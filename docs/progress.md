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
  - (A.19) の `Kb` と局所 amplitude root solver
  - (A.16) terminal `Qs` representation
  - (A.17) の `h^8`, `h^6` exact factors
  - CIで回帰し、診断CSV/日本語レポートを生成
- Appendix A.2 → A.3 の基準schedule積分
  - temporary reference power law → first transition → axial reduction → (A.9) hold を正規化ODEで接続
  - `m=M/(XE)`, `j=J/(XHE)`, `s=S/(XE^2)`, `rI=I/(XH)` を直接発展させる実装
  - raw `X,M,J,S` を生成せず巨大な対数半径を扱える
  - scheduleから得た `M,J` discrepancy を (A.15) two-bump closure に接続

### 数式監査・数値表現で検出した修正

- 類似座標の消去形は `q - z^2 q^(2h) = tau` と確認し、初期転記を修正した。
- Corollary A.3 のテストで power weights が重複する退化 exponent を避けるよう修正した。
- Appendix A.2 の巨大パラメータは `log P*`, `log h` を基本表現にした。
- Appendix A.3 の pulse moments は first bump center で正規化し、`exp(O(1/lambda))` を直接生成しない。
- A.11 の relative-E correction は局所 moment jump と global slope evolution を二重計上しないよう分離した。
- **重要な訂正:** (A.9) 後に列挙される4区間は A.3 closure のために使う補正ではない。Proposition A.4 はこれらを unaltered のまま残して (A.8) を成立させ、その後の cone realization / heat compensation / higher-order background / phase-averaged mean corrections が使用する。
- `lambda=0.05` を用いた初期 global test では `S/(XE^2)` が巨大化して A.19 bracket が消えた。これは reserved patch 不足ではなく、`Tw=60 log(1/lambda)` を含む ordered asymptotic construction に対して lambda が十分小さくなかったため。A.14 で使う `exp(-2 lambda Tw)=exp(-120 lambda log(1/lambda))` は lambda -> 0 で 1 に近づくが、lambda=0.05 では約 1.6e-8 まで落ちる。

### 成果物パイプライン

main 更新時に軽量基準計算を実行し、Taylor–Green検証、3D時系列VTK、可視化、レポートを GitHub に残す。Appendix A.3 は `results/reference/appendix_a3_closure.csv` と `reports/appendix-a3.md` を自動生成する。

## 現在地: Appendix A.3 global closure の漸近パラメータ検証

有限次元の A.3 closure mechanism と A.2 基準scheduleの正規化積分は実装済み。次に、原典の parameter choice order `(Md,Td,P*) -> lambda sufficiently small -> h sufficiently small` を数値テストにも反映し、full schedule から直接 (A.19) の root を回収できるか検証する。

A.9 の4つの reserved interval はこの段階では変更しない。原典どおり、A.3 / Proposition A.4 の成立後に後工程用の untouched patches として保持する。

次の順序:

1. **Appendix A.1 moment correction primitives — 実装済み**
2. **式 (4.15) cumulative moments — 実装済み**
3. **Corollary A.3 five-bump exact quadratic moment map — 実装済み**
4. **Appendix A.2 radial primitives / (A.10)–(A.13) — 実装済み**
5. **Appendix A.3 finite-dimensional closure mechanism — 実装済み**
6. **A.2 基準scheduleの正規化ODE — 実装済み**
7. **A.9 reserved correction intervals — 位置のみ実装済み、A.3では意図的に未変更**
8. ordered asymptotic regime で corrected A.14 / global A.19 bracket を数値確認
9. A.3 の `Amp(eta)` を global schedule から解き、(A.8) を全区間で確認
10. Appendix A.4 pressure datum `Pi0(eta)` / (A.21)–(A.23)
11. Lemma A.6 / Proposition A.7: exact heat exterior replacement と3-moment compensation
12. Appendix B analytic axis profile と moment matching
13. Appendix C admissible stress cone realization
14. Appendix C.2 後に first reserved patch で five moments を復元
15. higher-order / mean corrections で third/fourth reserved patches を使用
16. Theorem 4.6 profile assembly

任意の surrogate profile を論文構成そのものとして扱わない。

## その後

1. Sections 6–7 の oscillatory pulse
2. Propositions 7.5–7.6 の stress realization
3. 残差改善
4. Section 10 の localization と smooth forcing
5. 必要に応じて 64^3 → 128^3 → 256^3 解像度研究
