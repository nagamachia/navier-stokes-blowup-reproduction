# 解像度研究プロトコル

この文書は `docs/roadmap.md` の Step 4 を再現可能な数値手順として定義する。blow-up プロファイルを入れる前に、まず Taylor–Green 基準問題で数値基盤を検証する。

## 目的

解像度研究で答えたいのは「観測された差が空間解像度不足によるものか、それとも実装誤差・時間刻み誤差で説明できるか」である。

単に大きな計算が可能だからという理由で格子を増やさない。標準系列は次の通り。

1. `64^3`: デバッグと最初の収束基準
2. `128^3`: 64^3 が安定し、主要誤差が空間解像度由来と考えられる場合のみ
3. `256^3`: 64^3 と 128^3 の差から、さらに解像する意味が確認できた場合のみ

直接 DFT の小型テストは Fourier 規約の監査用参照実装として残し、実用的な格子では FFTW 版を使う。

## ビルド

Debian/Ubuntu:

```bash
sudo apt-get install libfftw3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

解像度実験を解釈する前に、CTest がすべて通っていることを確認する。

## FFTW Taylor–Green 実行系

```bash
./build/fftw_taylor_green <N> <final_time> <dt> <viscosity> <output.csv>
```

64^3 の最初の基準例:

```bash
./build/fftw_taylor_green 64 0.1 0.001 0.1 taylor_green_64.csv
```

128^3 の例:

```bash
./build/fftw_taylor_green 128 0.1 0.0005 0.1 taylor_green_128.csv
```

ただし、この2本だけを比較して空間収束を主張してはいけない。時間刻みも変わっているため、まず各 N で `dt` を半分にした計算を行い、時間刻み誤差が空間解像度差より十分小さいことを確認する。

## 現在の数値定式化

FFTW 実行系は以下を使用する。

- 各方向 `2*pi` の周期ボックス
- complex-to-complex FFTW
- forward 変換を総格子点数で正規化
- 整数 Fourier 波数によるスペクトル微分
- 非ゼロ波数に対する非圧縮射影 `P_k = I - k k^T / |k|^2`
- 成分ごとの 2/3 dealiasing
- 回転形式非線形項 `u x omega` の後に射影
- 4次 Runge–Kutta
- 粘性項 `-nu |k|^2 u_hat`

現在の Taylor–Green 初期条件は 2D Taylor–Green vortex を 3D ボックスへ埋め込んだもの。この場では射影後非線形項が消え、速度振幅は `exp(-2 nu t)` に従う。このため end-to-end 回帰には向くが、乱流的な解像度試験としてはまだ簡単である。

## CSV 診断量

`fftw_taylor_green` は各時間ステップについて次を出力する。

- `step`: ステップ番号
- `time`: 時刻
- `max_velocity`: 最大速度ノルム
- `max_vorticity`: 最大渦度ノルム
- `kinetic_energy`: `0.5 * sum |u_hat|^2`
- `enstrophy`: `0.5 * sum |omega_hat|^2`
- `viscous_dissipation`: `2 * nu * enstrophy`
- `divergence_l2`: `k dot u_hat` のスペクトル L2
- `pde_residual_l2`: `(u_n-u_{n-1})/dt - RHS(u_n)` の離散整合性残差
- `projected_nonlinear_l2`: 射影・dealiasing 後の回転形式非線形項 L2
- `dt`: 時間刻み
- `cfl`: `max_velocity * dt / dx`

blow-up 構成に固有の項ごとのキャンセルは、論文の exact profile と forcing を実装してから追加する。

## N を上げる前の受入条件

各解像度で少なくとも1回 `dt` 半減計算を行う。次の条件を満たしてから N を増やす。

1. CTest がすべて成功
2. `divergence_l2` が速度スペクトルに対して浮動小数点ノイズ近傍
3. CFL が明確に安定域内
4. `dt` 半減による診断量の変化が N 変更による差より十分小さい
5. 正規化・aliasing の明白なバグで診断量が増えていない
6. より大きな N を必要とする量に系統的な解像度依存が見える

Taylor–Green では最終速度を `exp(-2 nu t)` 倍した解析解とも比較する。内蔵回帰が許容誤差を超える場合、実行ファイルは非ゼロ終了する。

## GitHub 上の自動基準計算

CI では重い 64^3 計算を毎回走らせず、N=16 の軽量基準計算を実行する。目的はコード変更後にも診断量の形と可視化パイプラインが壊れていないことを確認すること。

生成物:

- `results/reference/taylor_green_n16.csv`
- `results/figures/taylor_green_energy.svg`
- `results/figures/taylor_green_errors.svg`
- `reports/latest.md`

64^3 以上の本格実験は、収束確認の段階で個人 PC または明示的な実験 workflow に分離する。

## Step 4 でまだ証明できないこと

Taylor–Green の収束は数値基盤の妥当性を示すだけで、有限時間特異性構成の再現ではない。論文構成を挿入する前に、`docs/theory.md` に記録した exact formula を順番に実装し、類似座標・微分作用素・leading profile・stress・oscillatory correction・forcing の各段階を個別に検証する。
