# Android / Web 可視化

## 目的

ParaView 用の VTK 成果物を GitHub 上に残したまま、Android の Chrome から Taylor–Green vortex の速度場を対話的に確認できるようにする。

## 構成

- ParaView 用 VTK: `results/reference/*.vtk`
- Web Viewer: `web/index.html`, `web/viewer.js`
- 公開: GitHub Pages (`.github/workflows/pages.yml`)
- 時系列の軽量確認: `results/figures/taylor_green_3d.gif`

Web Viewer は legacy ASCII VTK の `STRUCTURED_POINTS` と `VECTORS` を読み込み、速度ベクトルを WebGL で描画する。スマートフォンでは全格子点を描画すると重くなるため、既定では4点ごとに間引く。間引き幅とベクトル倍率は画面上で変更できる。

## Android での操作

GitHub Pages の URL を Chrome で開く。1本指ドラッグで回転し、ピンチで拡大・縮小する。

## GitHub Pages の初回設定

リポジトリの Settings → Pages → Build and deployment で Source を `GitHub Actions` にする。以降は `main` 更新時に `pages` workflow が Viewer と VTK/GIF を配置する。

リポジトリが private の場合、GitHub プランや組織ポリシーによって Pages の利用可否・公開範囲が異なる。Pages が利用できない場合でも VTK 自体は `results/reference/` に残り、ParaView から利用できる。

## 対応データ

現在の Viewer は以下を表示する。

- `taylor_green_numerical_n32.vtk`
- `taylor_green_reference_n16.vtk`

将来は 3D 時系列 VTK、渦度、Q-criterion、等値面を追加する。
