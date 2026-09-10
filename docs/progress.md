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
  - `q - z^2 q^(2h) = tau` の数値解法を実装
  - 既知の `(q,eta)` から `(tau,z)` を作り逆変換する回帰テストを実装
- Lemma 4.1, 式 (4.2): `T_b`, `Z_b`
  - 作用素を実装
  - manufactured profile を使い、物理座標の中心有限差分 `partial_t`, `partial_z` と比較
- 式 (4.6)–(4.7): `A_X`, `V0`, `Pi`
  - composite Simpson 積分による radial average
  - `U`, `partial_eta U` から非圧縮条件の `V0` を再構成
  - 軸正則形 `E=sqrt(2X)F` を使った pressure integral
  - polynomial / exponential manufactured profile で回帰
- 式 (4.8)–(4.11): leading stress
  - `docs/leading-stress.md` にPDFから式を抽出
  - `H,F,l,W,H_c,S_q,S_n,Q_s,N_s,p_s,a,b_s,T_0` を汎用 profile callback で実装
  - `Q_s/N_s` のODE恒等式、軸上極限、`T_0` の軸正則性を manufactured profile で回帰
- 式 (4.4)–(4.5): Cartesian 軸正則性
  - `E=sqrt(2X)F`, `V0=Xv0` を使った Cartesian leading field を実装
  - `r -> 0` で transverse velocity が `O(r)`、円柱成分を再構成できることを回帰
- 式 (4.15): five cumulative radial moments `M,I,J,S,Cp`
  - regular-axis 変数 `E=sqrt(2X)F` で正則に積分する実装を追加
  - polynomial manufactured profile で5成分を解析値と照合
- Appendix A.1, Lemma A.1: finitely many moment adjustments
  - separated radial bumps と異なる power weights から moment matrix `B_ij` を構成
  - partial-pivot Gaussian elimination と determinant 診断を実装
  - 非零 determinant と係数回収を回帰
- Appendix A.1, Lemma A.2: quadratic moment correction
  - `F_eta(c)=B(eta)c+Q_eta(c,c)` の zero-start fixed-point branch を実装
  - small quadratic manufactured system と `Q=0` の線形極限を回帰

### 数式監査で検出した修正

初期転記では消去形を `q - z^2/q^(2h) = tau` としていたが、式 (4.1) の

`z = q^(1/2-h) eta`

`tau = q(1-eta^2)`

からの代数消去、および論文 Section 4.1 の本文を再確認し、正しくは

`q - z^2 q^(2h) = tau`

であることを確認した。`T_b`, `Z_b` の有限差分照合がこの不整合を検出したため、座標ソルバと文書を修正した。

### 成果物パイプライン

main 更新時に軽量基準計算を実行し、Taylor–Green検証、3D時系列VTK、可視化、レポートを GitHub に残す。Web Viewer は Android から速度・渦度・Q-criterion・等値面・時間発展を確認できる。

## 現在地: Appendix A の constructive profile

Theorem 4.6 と Appendices A/B/C の constructive profile `E,U,Pi` の数値化に着手した。

現在、Appendix A.1 と Section 4.2 の共通 matching 基盤まで実装済み。次は Corollary A.3 を直接数値化し、2つの `U` bumps と3つの `E` bumps から5 momentsを同時に修正できることを確認する。

順序:

1. **Appendix A.1 moment correction primitives — 実装済み**
2. **式 (4.15) cumulative moments — 実装済み**
3. Corollary A.3 の 5-bump exact quadratic moment map — 次の対象
4. Appendix A.2 の staged outer profile (`U0,E0`) と parameter schedule
5. pressure datum `Pi0(eta)` と exterior power law
6. Appendix A.3–A.8 の moment closure / exact heat exterior replacement
7. Appendix B の analytic axis profile と moment matching
8. Appendix C の admissible stress cone realization
9. Theorem 4.6 profile assembly

任意の surrogate profile を OpenAI 構成として扱わない。

## その後

1. Sections 6–7 の oscillatory pulse
2. Propositions 7.5–7.6 の stress realization
3. 残差改善
4. Section 10 の localization と smooth forcing
5. 必要に応じて 64^3 → 128^3 → 256^3 解像度研究
