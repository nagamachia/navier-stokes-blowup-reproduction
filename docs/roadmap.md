# 実験ロードマップ

## 基本原則

計算資源を増やす前に、数学と実装の検証に時間を使う。前段階で観測された誤差が「実装ミス」ではなく「解像度不足」に起因すると判断できるまでは、格子点数を増やさない。

## Step 1 — 監査可能な理論抽出

成果物: `docs/theory.md`

後で実装する量について、必ず以下を記録する。

- 論文中の記号
- 節・式番号
- 次元付きか無次元か
- 特異時刻へ近づくときの漸近スケーリング
- 厳密式、漸近式、数値簡略化のどれか

推測したプロファイルを OpenAI 論文の構成として扱わない。

## Step 2 — 最小コア実験

成果物: `src/core_scaling.cpp`

最初は FFTW も Navier–Stokes ソルバも使わず、抽出した自己相似スケーリングだけを小さな計算で確認する。これにより、構成理解の誤りと PDE ソルバの誤りを分離する。

目標計算時間: 一般的な CPU で数秒〜数分。

## Step 3 — ソルバ検証

成果物: 周期擬スペクトルソルバと回帰テスト。

順序:

1. Fourier 微分
2. 発散ゼロ射影
3. 2/3 dealiasing
4. 既知 Fourier モードの粘性減衰
5. Taylor–Green vortex 回帰

この5項目が通るまで、blow-up 構成をソルバへ入れない。

現状、この Step は完了している。

## Step 4 — 解像度研究

標準系列:

- 64^3: デバッグと最初の基準
- 128^3: 通常の比較実験
- 256^3: 収束データが必要性を示した場合のみ

最低限追跡する量:

- 最大速度
- 最大渦度
- 運動エネルギー
- エンストロフィー
- 粘性散逸
- 発散誤差
- PDE 残差と主要項ごとのキャンセル
- 時間刻み
- CFL 数

現在は FFTW 版実行系と自動可視化が整っており、Taylor–Green 基準計算の自動生成を行っている。論文構成そのものの小規模回帰を先に進めるため、大規模 DNS はまだ主作業にしない。

## Step 4.5 — 論文構成の数値部品

大規模 DNS の前に、論文の構成そのものを小さく検証する。詳細な実装境界は `docs/progress.md` を正とする。

順序:

1. 式 (3.2)/(4.1) の類似座標 `q, eta, X` — **実装済み**
2. Lemma 4.1, 式 (4.2) の微分作用素 `T_b`, `Z_b` — **実装済み・有限差分照合あり**
3. 式 (4.6)–(4.7) の `A_X`, `V0`, `Pi` — **実装済み**
4. 軸上正則性 (4.4)–(4.5) — **実装済み**
5. 式 (4.8)–(4.11) の leading stress — **実装済み**
6. 式 (4.15) と Appendix A.1 の moment correction primitives — **実装済み**
7. Appendix A.2–A.3 の staged outer profile / finite-dimensional closure — **実装・CI回帰済み**
8. Appendix A.4 / Lemma A.5 pressure datum (A.21)–(A.23) — **実装・CI回帰済み**
9. Appendix A.6 / Lemma A.6 exact heat factor (A.32)–(A.39) — **実装・CI回帰済み**
10. Proposition A.7 second reserved patch の `Cp,S,I` 3-moment compensation core — **実装・CI回帰済み**
11. Proposition A.7 の full schedule coupling — **進行中**
    - terminal heat replacement の normalized discrepancy → second patch bridge — **実装・runtime回帰済み**
    - A.2/A.3 full schedule から `X_patch,e_patch,X_tail,e_tail` の相対スケール抽出 — **log-spaceで実装・回帰済み**
    - ordered regime の heat discrepancy を signed-log で評価 — **実装・回帰済み**
    - ordered regime の3 target の log spread と dominant linear correction を監査 — **実装・CI回帰済み**
    - quadratic term は dominant target より小さいが、最小 target より大きくなり得る精度階層を検出 — **実装・CI回帰済み**
    - exact heat factor と leading signed-log 式の相対誤差 bound — **実装・CI回帰済み**
    - 階層的 nonlinear correction / contraction を Lemma A.2 と対応する形で監査し、global `Cp,S,I` と axis pressure datum restoration の再現範囲を確定 — **次の対象**
12. Appendix B analytic axis profile と moment matching
13. Appendix C admissible stress cone realization
14. Sections 6–7 の oscillatory pulse と平均 stress
15. residual improvement と final forcing

A.9 の4つの reserved intervals は用途を混同しない。現在使用しているのは heat compensation 用の **second reserved patch** だけで、他は後段の construction 用に保持する。

## Step 5 — 動的リスケーリング

物理空間のコアが解像できなくなった場合は、単純に巨大 DNS へ進む前に動的リスケーリングを試す。縮小する構造を計算座標上で O(1) に保つことが目的。

## 計算資源ポリシー

標準予算は既存の個人 PC とオープンソースソフトウェア。GPU、クラウド HPC、MPI は、計測によって具体的なボトルネックが確認された後に導入する。

## 成果物ポリシー

途中成果も GitHub 上で追跡可能にする。

- CSV: `results/reference/`
- 可視化 SVG/GIF: `results/figures/`
- 自動レポート: `reports/`

Appendix A.6/A.7 の軽量診断は `results/reference/appendix_a6_a7_heat.csv` に保存し、heat ODE residual、3×3 compensation Jacobian、moment recovery error に加えて、ordered regime の signed-log target、target log spread、quadratic-to-dominant ratio、leading-asymptotic relative-error bound を追跡する。

main 更新時に GitHub Actions が軽量な基準計算を実行し、結果を自動コミットする。

## 科学的注意

数値計算は特異時刻直前の挙動やスケーリング則との整合性を示す証拠にはなるが、それだけで有限時間特異性を数学的に証明することはできない。逆に、低解像度で特異的挙動が見えないことも、解析構成が誤りである証拠にはならない。
