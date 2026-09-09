# Navier–Stokes Blow-up Reproduction

OpenAI が発表した 3 次元非圧縮 Navier–Stokes 方程式の有限時間特異性構成を、個人計算機で追跡可能な範囲から段階的に数値検証するプロジェクトです。

> [!IMPORTANT]
> このリポジトリは独立した数値実験です。数学的証明、形式検証、あるいは発表結果の独立な厳密検証を主張するものではありません。目的は、特異時刻直前のスケーリング、座標変換、残差キャンセル構造を数値的に追跡することです。

## 目標

1. 論文から数値実装に必要な定義・式・仮定を式番号付きで抽出する。
2. 完全な PDE ソルバの前に、自己相似コアのスケーリングと類似座標を検証する。
3. 小型の 3D 周期擬スペクトル Navier–Stokes ソルバを構築・検証する。
4. 64^3、128^3、必要なら 256^3 の解像度比較を行う。
5. 最大速度、最大渦度、運動エネルギー、エンストロフィー、散逸、発散誤差、PDE 残差、CFL を追跡する。
6. 直接 DNS が非効率になった場合は動的リスケーリングを検討する。

## 設計方針

- CPU 優先
- C++20 / FP64
- FFTW を使用
- まず単一ノード・個人 PC で検証
- GPU / MPI / クラウド HPC は、具体的なボトルネックが確認されてから導入
- 解像度を上げる前に、低解像度で実装誤差・時間刻み誤差を潰す

## 現在の論文構成実装

### 類似座標 — 式 (3.2)/(4.1)

`tau = q (1 - eta^2)`

`z = q^(1/2-h) eta`

`X = r^2 / (2q)`

を扱い、`include/similarity_coordinates.hpp` で

`q - z^2 q^(2h) = tau`

を数値的に解きます。逆変換回帰では `q,eta,X` を元の恒等式へ戻して検証します。

### 微分作用素 — Lemma 4.1 / 式 (4.2)

`include/similarity_operators.hpp` に `T_b`, `Z_b` を実装し、manufactured profile を使って物理座標の有限差分 `partial_t`, `partial_z` と照合します。

### Leading profile 関係 — 式 (4.6)–(4.7)

`include/leading_profile_relations.hpp` に

- radial average `A_X`
- `U`, `partial_eta U` からの `V0`
- 軸正則形 `E=sqrt(2X)F` を使った `Pi`

を実装し、解析値を持つ manufactured profile と比較します。

### Leading stress — 式 (4.8)–(4.11)

`docs/leading-stress.md` にPDFから式を抽出し、`include/leading_stress.hpp` に

`H,F,l,W,H_c,S_q,S_n,Q_s,N_s,p_s,a,b_s,T_0`

を実装しています。`Q_s/N_s` のODE恒等式、軸上極限、`T_0` の軸正則性を回帰テストしています。

### Cartesian 軸正則性 — 式 (4.4)–(4.5)

`include/leading_field.hpp` に `E=sqrt(2X)F`, `V0=Xv0` を使った Cartesian leading field を実装し、`r -> 0` で transverse velocity が `O(r)` になることを検証します。

## 擬スペクトルソルバ検証

CTest では Fourier 微分、発散ゼロ射影、2/3 dealiasing、粘性減衰、Taylor–Green、FFTW smoke test と上記論文由来の数値部品をまとめて回帰検証します。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Debian/Ubuntu では事前に

```bash
sudo apt-get install libfftw3-dev
```

が必要です。

## 実験結果を見る

main 更新時に GitHub Actions が軽量な基準計算を自動実行し、結果をリポジトリへコミットします。

- 生データ: `results/reference/`
- 可視化: `results/figures/`
- 最新レポート: [reports/latest.md](reports/latest.md)

Taylor–Green は `z=0` 断面について、流れそのものも可視化します。

- `results/figures/taylor_green_velocity_field.svg` — 速度の大きさ + ベクトル
- `results/figures/taylor_green_streamlines.svg` — 流線
- `results/figures/taylor_green_vorticity.svg` — 渦度 `omega_z`
- `results/figures/taylor_green_energy.svg` — エネルギー・エンストロフィー
- `results/figures/taylor_green_errors.svg` — 数値誤差診断
- `results/figures/similarity_coordinate_error.svg` — 類似座標誤差

最新レポートにも速度場・流線・渦度図を埋め込みます。

## FFTW 解像度実験

64^3 の Taylor–Green 基準計算例:

```bash
./build/fftw_taylor_green 64 0.1 0.001 0.1 taylor_green_64.csv
```

CSV には最大速度、最大渦度、運動エネルギー、エンストロフィー、粘性散逸、発散 L2、離散 PDE 残差、射影後非線形項 L2、時間刻み、CFL を出力します。

64^3 → 128^3 → 256^3 の昇格条件は `docs/resolution-study.md` に記載しています。

## ドキュメント

- `docs/theory.md`: 論文から抽出した数式・式番号・実装境界
- `docs/leading-stress.md`: 式 (4.8)–(4.11) の stress 定義
- `docs/roadmap.md`: 開発順序
- `docs/resolution-study.md`: 解像度実験手順
- `docs/progress.md`: 実装済み範囲と次の対象
- `reports/latest.md`: CI が生成する最新数値レポート

説明文書と自動レポートは原則として日本語で管理します。

## 現在の進捗

- **Step 1:** 類似座標、微分作用素、leading field、非圧縮条件、圧力、leading stress、active annulus、最終 forcing 構成まで式番号付きで抽出を継続中。
- **Step 2:** コアスケーリング実験済み。
- **Step 3:** 擬スペクトルソルバの5段階検証を完了。
- **Step 4 基盤:** FFTW 版 3D 実行系、CSV診断、SVG可視化、自動レポート生成を実装済み。
- **論文構成の汎用数値部品:** 式 (3.2)/(4.1)、(4.2)、(4.4)–(4.11) まで実装・回帰中。
- **次の大きな対象:** Theorem 4.6 と Appendices A/B の constructive profile `E,U,Pi`。ここからは実際の profile construction を数値化し、任意の surrogate を OpenAI 構成として扱わない。

## 原典

- OpenAI 発表: https://openai.com/index/navier-stokes-solution/
- 論文: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- 形式化リポジトリ: https://github.com/openai/NavierStokesAndEuler
