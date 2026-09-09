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
  - `q - z^2/q^(2h) = tau` の数値解法を実装
  - 既知の `(q,eta)` から `(tau,z)` を作り逆変換する回帰テストを実装
- Lemma 4.1, 式 (4.2): `T_b`, `Z_b`
  - 作用素を実装
  - manufactured profile を使い、物理座標の中心有限差分 `partial_t`, `partial_z` と比較する回帰テストを実装

### 成果物パイプライン

main 更新時に軽量基準計算を実行し、次を GitHub に残す。

- `results/reference/taylor_green_n16.csv`
- `results/reference/similarity_coordinates.csv`
- `results/figures/taylor_green_energy.svg`
- `results/figures/taylor_green_errors.svg`
- `results/figures/similarity_coordinate_error.svg`
- `reports/latest.md`

## 次の実装対象

式 (4.6)–(4.7) の

- radial average `A_X`
- 非圧縮条件から決まる `V0`
- pressure relation `partial_X Pi = E^2/(2X)`

を、まず manufactured profile で検証可能な数値部品として実装する。

ただし exact `E,U` は単純な閉形式ではないため、Theorem 4.6 と Appendices A/B の constructive profile を無視して適当な代替関数を OpenAI 構成として扱わない。

## その後

1. 軸上正則性 (4.4)–(4.5)
2. 式 (4.8)–(4.11) の leading stress の詳細抽出
3. active annulus
4. Sections 6–7 の oscillatory pulse
5. Propositions 7.5–7.6 の stress realization
6. 残差改善
7. Section 10 の localization と smooth forcing
8. 必要に応じて 64^3 → 128^3 → 256^3 解像度研究
