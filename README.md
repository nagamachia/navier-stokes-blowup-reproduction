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

## 現在の構成

### 論文由来の類似座標

OpenAI 論文の式 (3.2)/(4.1) に対応して

`tau = q (1 - eta^2)`

`z = q^(1/2-h) eta`

`X = r^2 / (2q)`

を扱います。`include/similarity_coordinates.hpp` には、与えられた `(z, tau)` から

`q - z^2 q^(2h) = tau`

を数値的に解く実装があります。`tests/similarity_coordinates.cpp` で元の恒等式へ戻ることを回帰テストしています。

### Lemma 4.1 の微分作用素

`include/similarity_operators.hpp` に式 (4.2) の `T_b`, `Z_b` を実装しています。`tests/similarity_operators.cpp` では manufactured profile を使い、物理座標の有限差分 `partial_t`, `partial_z` と照合します。

### 式 (4.6)–(4.7) の leading profile 関係

`include/leading_profile_relations.hpp` に以下の汎用数値部品を実装しています。

- radial average `A_X`
- `U`, `partial_eta U` からの `V0`
- 軸正則形 `E=sqrt(2X)F` を使った pressure integral `Pi`

`tests/leading_profile_relations.cpp` では polynomial / exponential manufactured profile に対して解析値と比較します。

### Leading stress

式 (4.8)–(4.11) の `H,F,l,W,H_c,S_q,S_n,Q_s,N_s,T_0` は、論文PDFから直接抽出して `docs/leading-stress.md` に整理済みです。次はこの部分の数値実装です。

### 擬スペクトルソルバ検証

CTest では Fourier 微分、発散ゼロ射影、2/3 dealiasing、既知 Fourier モードの粘性減衰、Taylor–Green vortex、FFTW smoke test、類似座標、類似微分作用素、leading-profile relation を回帰検証します。

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

GitHub Actions は main 更新時に小さな基準計算を自動実行し、結果をリポジトリへ保存します。

- 生データ: `results/reference/`
- 図: `results/figures/`
- 最新レポート: `reports/latest.md`

GitHub 上で直接確認できる図:

- `results/figures/taylor_green_energy.svg`
- `results/figures/taylor_green_errors.svg`
- `results/figures/similarity_coordinate_error.svg`

最新のまとめは [reports/latest.md](reports/latest.md) を参照してください。

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

- **Step 1:** 類似座標、微分作用素、leading field、非圧縮条件、圧力、leading stress、active annulus、最終 forcing 構成まで式番号付きで抽出中。
- **Step 2:** コアスケーリング実験済み。
- **Step 3:** Fourier 微分、射影、2/3 dealiasing、粘性減衰、Taylor–Green 回帰が CI で通過。
- **Step 4 基盤:** FFTW 版 3D 実行系、診断量、CSV、SVG可視化、自動レポート生成を実装済み。
- **論文構成の実装:** 式 (3.2)/(4.1)、Lemma 4.1、式 (4.6)–(4.7) の汎用数値部品まで実装済み。次は式 (4.8)–(4.11) の leading stress を実装する。

## 原典

- OpenAI 発表: https://openai.com/index/navier-stokes-solution/
- 論文: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- 形式化リポジトリ: https://github.com/openai/NavierStokesAndEuler
