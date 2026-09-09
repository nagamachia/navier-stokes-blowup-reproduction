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

`include/similarity_operators.hpp` に式 (4.2) の

`T_b f`

`Z_b f`

を実装しています。`tests/similarity_operators.cpp` では、任意の manufactured profile を使って、物理座標で直接計算した有限差分 `partial_t`, `partial_z` と照合します。

### 擬スペクトルソルバ検証

CTest では以下を検証しています。

1. Fourier 微分
2. 発散ゼロ射影
3. 2/3 dealiasing
4. 既知 Fourier モードの粘性減衰
5. Taylor–Green vortex 回帰
6. FFTW 版 Taylor–Green smoke test
7. 類似座標の逆変換回帰
8. Lemma 4.1 の `T_b`, `Z_b` 有限差分回帰

実行:

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

特に以下の図を GitHub 上で直接確認できます。

- `results/figures/taylor_green_energy.svg`
- `results/figures/taylor_green_errors.svg`
- `results/figures/similarity_coordinate_error.svg`

最新のまとめは [reports/latest.md](reports/latest.md) を参照してください。

## FFTW 解像度実験

例として 64^3 の Taylor–Green 基準計算は

```bash
./build/fftw_taylor_green 64 0.1 0.001 0.1 taylor_green_64.csv
```

で実行できます。

CSV には以下を出力します。

- 最大速度
- 最大渦度
- 運動エネルギー
- エンストロフィー
- 粘性散逸
- 発散 L2 誤差
- 離散 PDE 残差
- 射影後非線形項 L2
- 時間刻み
- CFL

64^3 → 128^3 → 256^3 の昇格条件は `docs/resolution-study.md` に記載しています。

## ドキュメント

- `docs/theory.md`: 論文から抽出した数式・式番号・実装境界
- `docs/roadmap.md`: 開発順序
- `docs/resolution-study.md`: 解像度実験手順
- `docs/progress.md`: 実装済み範囲と次の対象
- `reports/latest.md`: CI が生成する最新数値レポート

以後、説明文書と自動レポートは原則として日本語で管理します。

## 現在の進捗

- **Step 1:** 論文から類似座標、微分作用素、leading field、非圧縮条件、圧力、active annulus、最終 forcing 構成まで式番号付きで抽出済み。
- **Step 2:** コアスケーリング実験済み。
- **Step 3:** Fourier 微分、射影、2/3 dealiasing、粘性減衰、Taylor–Green 回帰が CI で通過。
- **Step 4 基盤:** FFTW 版 3D 実行系、診断量、CSV 出力、可視化・レポート自動生成まで実装済み。
- **論文構成の実装:** 式 (3.2)/(4.1) の類似座標と Lemma 4.1 の `T_b`, `Z_b` を実装済み。次は式 (4.6)–(4.7) の radial average `A_X`、`V0`、`Pi` を実装する。

## 原典

- OpenAI 発表: https://openai.com/index/navier-stokes-solution/
- 論文: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- 形式化リポジトリ: https://github.com/openai/NavierStokesAndEuler
