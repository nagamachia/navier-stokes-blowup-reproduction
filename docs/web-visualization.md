# Android / Web 可視化

## 目的

ParaView 用の VTK 成果物を GitHub 上に残したまま、Android の Chrome から Taylor–Green vortex の速度場と渦構造を対話的に確認できるようにする。

## 構成

- ParaView 用 VTK: `results/reference/*.vtk`
- Web Viewer: `web/index.html`, `web/viewer.js`
- 公開: GitHub Pages (`.github/workflows/pages.yml`)
- 時系列の軽量確認: `results/figures/taylor_green_3d.gif`

Web Viewer は legacy ASCII VTK の `STRUCTURED_POINTS` と `VECTORS` を読み込み、WebGL で描画する。スマートフォンでは全格子点を描画すると重くなるため、既定では4点ごとに間引く。

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

渦度と Q-criterion の「強い領域の表示割合」は画面から変更できる。例えば20%なら対象スカラーの上位20%を表示する。Q-criterion では正値だけを母集団にする。

この計算は可視化用の軽量診断であり、論文用の高精度微分を置き換えるものではない。研究用の定量値はソルバ側の Fourier 微分や保存済み診断量を優先する。

## Android での操作

GitHub Pages の URL を Chrome で開く。1本指ドラッグで回転し、ピンチで拡大・縮小する。

公開 URL:

`https://nagamachia.github.io/navier-stokes-blowup-reproduction/`

## GitHub Pages の初回設定

リポジトリの Settings → Pages → Build and deployment で Source を `GitHub Actions` にする。以降は `main` 更新時に `pages` workflow が Viewer と VTK/GIF を配置する。

## 対応データ

現在の Viewer は以下を表示する。

- `taylor_green_numerical_n32.vtk`
- `taylor_green_reference_n16.vtk`

## 次の拡張

- 3D 時系列 VTK と再生スライダー
- 渦度・Q-criterion の等値面
- 時刻ごとの `max |ω|`, energy, enstrophy と3D表示の同期
- 解像度系列 64^3 / 128^3 の比較表示
