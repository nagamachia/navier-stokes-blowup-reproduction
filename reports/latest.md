# 最新の自動計算レポート

このレポートは GitHub Actions で自動生成されます。対象コミット: `f849d294205d`

## 実施した計算

- FFTW 版 Taylor–Green 基準計算
- OpenAI 論文の式 (3.2)/(4.1) に対応する類似座標 `q, eta, X` の再構成テスト
- Lemma 4.1, 式 (4.2) の `T_b`, `Z_b` と物理座標有限差分の照合
- Fourier 微分・非圧縮射影・2/3 dealiasing・粘性減衰・Taylor–Green 回帰

## Taylor–Green 基準計算の要約

- 最終時刻: `0.1`
- 最終運動エネルギー: `0.240197359788`
- 最終エンストロフィー: `0.480394719576`
- 最大発散誤差 L2: `4.567e-32`
- 最大 CFL: `0.0127324`

![Taylor–Green energy](../results/figures/taylor_green_energy.svg)

![Taylor–Green errors](../results/figures/taylor_green_errors.svg)

## 類似座標ソルバの要約

論文の恒等式

`tau = q (1 - eta^2)`

`z = q^(1/2-h) eta`

を使って既知の `(q, eta)` から `(tau, z)` を作り、数値ソルバで `q` と `eta` を逆算しています。

- 最大 `q` 相対誤差: `8.882e-15`
- 最大 `eta` 絶対誤差: `3.442e-15`

![Similarity coordinate errors](../results/figures/similarity_coordinate_error.svg)

## 解釈

現在は、論文固有の座標変換と Lemma 4.1 の微分作用素まで実装しています。次の実装対象は式 (4.6)–(4.7) に基づく radial average `A_X`、`V0`、`Pi` の再構成です。exact profile `E,U` 自体は Theorem 4.6 と Appendices A/B の構成を数値化してから導入します。

CSV の生データは `results/reference/` に保存されています。
