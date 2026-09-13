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
11. Proposition A.7 full schedule coupling / nonlinear recovery audit — **数値再現完了・CI回帰済み**
    - terminal heat replacement の normalized discrepancy → second patch bridge
    - A.2/A.3 full schedule から `X_patch,e_patch,X_tail,e_tail` の相対スケールを log-space で抽出
    - ordered regime の heat discrepancy を signed-log で評価
    - exact heat factor と leading signed-log 式の相対誤差を上から評価
    - 3 target の巨大な log spread を明示し、通常倍精度での componentwise exact solve を禁止
    - `B^{-1}` と quadratic map `Q` の有限次元 operator bound から Lemma A.2 fixed-point map の contraction factor を log-space で評価
    - exact heat remainder を contraction bound に伝播し、exact heat discrepancy に対する small nonlinear recovery branch の存在・一意性を監査
    - `Cp` が axis pressure datum increment と同じ条件であることを end-to-end result に統合
12. Appendix B analytic axis profile と moment matching — **進行中**
    - B.1 の `U*`,`H*`,`W*`,`Z*` を A.4 の `Pi0,Pi0_eta` に接続 — **実装・CI回帰済み**
    - `H*` の一意な零点 `eta0 in (-1,0)`、`Z*(eta0)>0`、`-W*>2.8` — **数値監査済み**
    - B.2 の small-`|Z*|` set と `chi>0.99` separation — **実装・回帰済み**
    - B.3 の `phi*`, B.11 の比較関数 `f0` と `f0>0.265` — **実装・回帰済み**
    - B.13 の比較解 `Phi0=f0(Y chi)`, `u0=-YZ*/(2L)` と B.19 endpoint の二分岐機構 — **漸近監査・CI回帰済み**
    - `chi<=0.99` 領域は `sigma*` が非常に小さいため uniform eta grid が見落とし得ることを検出し、`H*(eta0)=0` の根を endpoint audit に明示的に含めるよう修正
    - B.5 の coefficient operator `J1,J2`、`T=(1/2)J2(chi·)`、有限次数 `(1+T)^-1` と `Phi0` coefficient series — **実装・CI回帰済み**
    - weighted coefficient l1 norm と `Jnu` / `T` の operator bound — **実装・CI回帰済み**
    - `L^{-1}(g + Lambda^{-1} Q(z,z))` 型の Banach fixed-point contraction audit scaffold — **実装・CI回帰中**
    - 原典 Proposition B.2 の nonlinear map を式番号どおりに coefficient-space へ落とし、上記 scaffold に実係数を接続 — **次の対象**
    - 後続 continuation・five-moment matching — **未着手**
13. Appendix C admissible stress cone realization
14. Sections 6–7 の oscillatory pulse と平均 stress
15. residual improvement と final forcing

A.9 の4つの reserved intervals は用途を混同しない。現在使用済みなのは heat compensation 用の **second reserved patch** だけで、他は後段の construction 用に保持する。

## Step 5 — 動的リスケーリング

物理空間のコアが解像できなくなった場合は、単純に巨大 DNS へ進む前に動的リスケーリングを試す。縮小する構造を計算座標上で O(1) に保つことが目的。

## 計算資源ポリシー

標準予算は既存の個人 PC とオープンソースソフトウェア。GPU、クラウド HPC、MPI は、計測によって具体的なボトルネックが確認された後に導入する。

## 成果物ポリシー

途中成果も GitHub 上で追跡可能にする。

- CSV: `results/reference/`
- 可視化 SVG/GIF: `results/figures/`
- 自動レポート: `reports/`

Appendix A.6/A.7 の軽量診断は `results/reference/appendix_a6_a7_heat.csv` に保存し、heat ODE residual、3×3 compensation Jacobian、moderate-regime moment recovery error、ordered-regime signed-log target、target log spread、quadratic scale、leading-asymptotic relative-error bound、nonlinear contraction bound、exact-target recovery certification を追跡する。

main 更新時に GitHub Actions が軽量な基準計算を実行し、結果を自動コミットする。

## 科学的注意

ここでの「A.7 数値再現完了」は、通常倍精度で全ての指数的に小さい係数を直接表示したという意味ではない。moderate regime では直接回帰し、ordered regime では log-space の厳密なスケール管理と有限次元 contraction bound により再現範囲を閉じている。これは Proposition A.7 の数学的証明そのものを新たに与えるものではない。

Appendix B でも、比較解・漸近式の回帰と full nonlinear analytic construction を区別する。B.13/B.19 の比較機構や B.5 の線形 coefficient operator、一般形の contraction scaffold が数値的に確認できても、それだけで Proposition B.2 の nonlinear coefficient-space contraction や後続 moment matching を完了扱いにはしない。原典の nonlinear map の係数を抽出して接続した時点で初めて B.2 の full audit とする。

数値計算は特異時刻直前の挙動やスケーリング則との整合性を示す証拠にはなるが、それだけで有限時間特異性を数学的に証明することはできない。逆に、低解像度で特異的挙動が見えないことも、解析構成が誤りである証拠にはならない。
