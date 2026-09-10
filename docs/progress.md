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
  - 5×5 Jacobian を数値構成
  - momentごとのスケール差を行正規化
  - 線形係数回収と非線形 fixed-point 回収を回帰
- Appendix A.2 staged outer profile の radial primitives
  - (A.5) smooth step `sigma`
  - (A.6) parameter hierarchy を `log P*`, `log h` で安全に表現
  - (A.7) reference inner power law
  - axial reduction / (A.9) intermediate power law / reserved patches
  - (A.10) profile interpolation
  - (A.12) terminal factor / slope
  - (A.13) `Qp`
- Appendix A.3 moment closure mechanism
  - axial pulse `phi_b`, `R0`
  - (A.15) の2-bump `M=J=0` closure
  - (A.11) の2-bump `I` / pressure increment closure
  - (A.19) の `Kb` と amplitude root solver
  - (A.16) terminal `Qs` representation
  - (A.17) の `h^8`, `h^6` exact factors
  - CIで回帰し、診断CSV/日本語レポートを生成

### 数式監査・数値表現で検出した修正

- 類似座標の消去形は `q - z^2 q^(2h) = tau` と確認し、初期転記を修正した。
- Corollary A.3 のテストで power weights が重複する退化 exponent を避けるよう修正した。
- Appendix A.2 の parameter hierarchy `Td=exp(Md)+10`, `P*>exp(Td)`, `h<exp(-Td)` は通常の `double` で `P*`,`h` を直接保持するとオーバーフロー/アンダーフローするため、実装では `log P*`, `log h` を基本表現にした。
- Appendix A.3 の pulse moments は first bump center で正規化し、`exp(O(1/lambda))` を直接生成しない実装にした。

### 成果物パイプライン

main 更新時に軽量基準計算を実行し、Taylor–Green検証、3D時系列VTK、可視化、レポートを GitHub に残す。Web Viewer は Android から速度・渦度・Q-criterion・等値面・時間発展を確認できる。

Appendix A.3 については `results/reference/appendix_a3_closure.csv` と `reports/appendix-a3.md` を自動生成する。

## 現在地: Appendix A.3 統合 → A.4 pressure datum

Appendix A.3 の有限次元 closure mechanism は実装・回帰済み。ただし、現状の amplitude solver は (A.19) の小さい remainder `Error(Amp,eta)` を callback として受け取る段階で、A.2 の全radial scheduleからそのremainderを直接計算する統合は未完。

次の順序:

1. **Appendix A.1 moment correction primitives — 実装済み**
2. **式 (4.15) cumulative moments — 実装済み**
3. **Corollary A.3 five-bump exact quadratic moment map — 実装済み**
4. **Appendix A.2 radial primitives / (A.10)–(A.13) — 実装済み**
5. **Appendix A.3 finite-dimensional closure mechanism — 実装済み**
6. A.2 全stageをglobal log-radial scheduleとして接続し、(A.14)/(A.19) remainderを直接評価
7. A.3 の `Amp(eta)`, pulse corrections, angular correctionsを実profileへ統合して (A.8) を全区間で確認
8. Appendix A.4 pressure datum `Pi0(eta)` / (A.21)–(A.23)
9. Lemma A.6 / Proposition A.7: exact heat exterior replacement と3-moment compensation
10. Appendix B analytic axis profile と moment matching
11. Appendix C admissible stress cone realization
12. Theorem 4.6 profile assembly

任意の surrogate profile を OpenAI 構成として扱わない。

## その後

1. Sections 6–7 の oscillatory pulse
2. Propositions 7.5–7.6 の stress realization
3. 残差改善
4. Section 10 の localization と smooth forcing
5. 必要に応じて 64^3 → 128^3 → 256^3 解像度研究
