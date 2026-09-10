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
  - temporary reference power law → first transition → axial reduction → 未補正の (A.9) hold を正規化ODEで接続
  - `m=M/(XE)`, `j=J/(XHE)`, `s=S/(XE^2)`, `rI=I/(XH)` を直接発展させる実装
  - raw `X,M,J,S` を生成せず巨大な対数半径を扱える
  - 未補正scheduleから得た `M,J` discrepancy を (A.15) two-bump closure に接続

### 数式監査・数値表現で検出した修正

- 類似座標の消去形は `q - z^2 q^(2h) = tau` と確認し、初期転記を修正した。
- Corollary A.3 のテストで power weights が重複する退化 exponent を避けるよう修正した。
- Appendix A.2 の巨大パラメータは `log P*`, `log h` を基本表現にした。
- Appendix A.3 の pulse moments は first bump center で正規化し、`exp(O(1/lambda))` を直接生成しない。
- A.11 の relative-E correction は局所 moment jump と global slope evolution を二重計上しないよう分離した。
- **重要:** 未補正の (A.9) hold をそのまま積分すると `S/(XE^2)` が大きな負値へ増幅し、`S(infinity)=0` の amplitude root は `[0.9,1.2]` に現れない。これは論文が (A.9) 内に予約した correction patches をまだ適用していないためであり、A.14 の完成状態として扱ってはいけない。

### 成果物パイプライン

main 更新時に軽量基準計算を実行し、Taylor–Green検証、3D時系列VTK、可視化、レポートを GitHub に残す。Appendix A.3 は `results/reference/appendix_a3_closure.csv` と `reports/appendix-a3.md` を自動生成する。

## 現在地: Appendix A.9 reserved corrections → A.3 global closure

有限次元の A.3 closure mechanism は実装済み。global schedule については基準stageの正規化積分まで実装したが、A.9 の4つの reserved correction patches をまだ profile に反映していない。

CIでこの欠落を実際に検出した。未補正scheduleでは `S(infinity)` が amplitude 0.9 と1.2の両方で約 `-9.83e8` となり root が消えるため、full A.3 closure の成功条件から一旦外し、この値を「予約patch未適用の診断」として扱う。

次の順序:

1. **Appendix A.1 moment correction primitives — 実装済み**
2. **式 (4.15) cumulative moments — 実装済み**
3. **Corollary A.3 five-bump exact quadratic moment map — 実装済み**
4. **Appendix A.2 radial primitives / (A.10)–(A.13) — 実装済み**
5. **Appendix A.3 finite-dimensional closure mechanism — 実装済み**
6. **A.2 基準scheduleの正規化ODE — 実装済み**
7. **A.9 reserved correction intervals — 位置のみ実装済み**
8. A.9 の reserved profile/heat/positivity/mean corrections を原典どおり適用し、corrected A.14 state を構成
9. corrected A.14 → pulse → A.10 → A.11 → exterior transition → terminal tail を接続
10. global `S(infinity)` の `[0.9,1.2]` root と `Amp(eta)` を直接確認し、(A.8) を全区間で監査
11. Appendix A.4 pressure datum `Pi0(eta)` / (A.21)–(A.23)
12. Lemma A.6 / Proposition A.7: exact heat exterior replacement と3-moment compensation
13. Appendix B analytic axis profile と moment matching
14. Appendix C admissible stress cone realization
15. Theorem 4.6 profile assembly

任意の surrogate profile を論文構成そのものとして扱わない。

## その後

1. Sections 6–7 の oscillatory pulse
2. Propositions 7.5–7.6 の stress realization
3. 残差改善
4. Section 10 の localization と smooth forcing
5. 必要に応じて 64^3 → 128^3 → 256^3 解像度研究
