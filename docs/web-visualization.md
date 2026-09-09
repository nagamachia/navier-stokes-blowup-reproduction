# Android / Web 可視化

## 目的

ParaView 用の VTK 成果物を GitHub 上に残したまま、Android の Chrome から Taylor–Green vortex の速度場・渦構造・時間発展を対話的に確認できるようにする。

## 構成

- 静的 ParaView 用 VTK: `results/reference/*.vtk`
- 3D 時系列 VTK: `results/reference/taylor_green_3d_series/*.vtk`
- 時系列 manifest: `results/reference/taylor_green_3d_series/manifest.json`
- Web Viewer: `web/index.html`, `web/viewer.js`
- 公開: GitHub Pages (`.github/workflows/pages.yml`)
- 軽量 GIF: `results/figures/taylor_green_3d_evolution.gif`

3D Taylor–Green ソルバは従来の中央 `xy/xz/yz` 断面 CSV に加えて、オプションで全3D格子を各保存時刻に出力する。CI では一時的な全格子 CSV を `scripts/generate_taylor_green_3d_series.py` で時刻別 legacy ASCII VTK に変換し、11フレームの manifest とともに GitHub に保存する。

現在の標準 Web 時系列は `N=24`, `t=0.0 ... 0.5`, 11フレーム。全格子CSVはリポジトリには残さず、時刻別VTKのみを追跡する。

## 表示モード

現在は次の3種類を切り替えられる。

1. 速度ベクトル `|u|`
   - 線分で速度方向と大きさを表示する。
   - 色は速度の大きさで変化する。
2. 渦度強度 `|ω|`
   - `ω = curl u` をブラウザ側で周期中心差分から計算する。
   - 強い領域を点群として表示する。
3. Q-criterion
   - 速度勾配 `A = grad u` から対称部分 `S` と反対称部分 `Omega` を計算し、`Q = 1/2 (||Omega||^2 - ||S||^2)` を評価する。
   - `Q > 0` のうち強い領域を点群表示する。

渦度と Q-criterion の「強い領域の表示割合」は画面から変更できる。Q-criterion では正値だけを母集団にする。

このブラウザ側微分は可視化用の軽量診断であり、研究用の定量値はソルバ側の Fourier 微分や保存済み診断量を優先する。

## 時間発展

データ欄で `3D Taylor–Green N=24 time series` を選ぶと、`manifest.json` を読み込み、時刻スライダーと再生ボタンが有効になる。

- スライダー: 任意フレームへ移動
- `再生`: 約0.65秒間隔で次フレームへ進む
- 最終フレーム後は先頭へ戻る
- カメラ位置は時間切り替え中に保持する
- 直近のVTKフレームを小さくキャッシュし、Androidでの再読込を抑える

## Android での操作

GitHub Pages の URL を Chrome で開く。1本指ドラッグで回転し、ピンチで拡大・縮小する。

公開 URL:

`https://nagamachia.github.io/navier-stokes-blowup-reproduction/`

## GitHub Pages

リポジトリの Settings → Pages → Build and deployment で Source を `GitHub Actions` にする。`main` 更新時に `pages` workflow が Viewer、静的 VTK、3D 時系列 VTK、GIF を配置する。

## 対応データ

- `taylor_green_3d_series/manifest.json` + `frame_000.vtk` ... `frame_010.vtk`
- `taylor_green_numerical_n32.vtk`
- `taylor_green_reference_n16.vtk`

## 次の拡張

- 渦度・Q-criterion の等値面
- 時刻ごとの `max |ω|`, energy, enstrophy と3D表示の同期
- 解像度系列 64^3 / 128^3 の比較表示
- VTK ASCII より軽量なバイナリ配信形式の検討
