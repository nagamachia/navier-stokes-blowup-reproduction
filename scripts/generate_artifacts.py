#!/usr/bin/env python3

import csv
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures"
REPORT_DIR = ROOT / "reports"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"CSV is empty: {path}")
    return rows


def f(row, key):
    try:
        return float(row[key])
    except (KeyError, ValueError):
        return math.nan


def finite(values):
    return [v for v in values if math.isfinite(v)]


def svg_line_plot(path, title, x_values, series, x_label, y_label, log_y=False):
    width, height = 920, 520
    left, right, top, bottom = 90, 30, 65, 75
    plot_w = width - left - right
    plot_h = height - top - bottom

    x_vals = finite(x_values)
    if not x_vals:
        raise RuntimeError("no finite x values")
    xmin, xmax = min(x_vals), max(x_vals)
    if xmin == xmax:
        xmax = xmin + 1.0

    transformed = []
    for _, values in series:
        for value in values:
            if not math.isfinite(value):
                continue
            transformed.append(
                math.log10(max(abs(value), 1e-18)) if log_y else value
            )
    if not transformed:
        transformed = [0.0, 1.0]
    ymin, ymax = min(transformed), max(transformed)
    if ymin == ymax:
        pad = max(1.0, abs(ymin) * 0.1)
        ymin -= pad
        ymax += pad
    else:
        pad = 0.08 * (ymax - ymin)
        ymin -= pad
        ymax += pad

    def sx(x):
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def sy(y):
        value = math.log10(max(abs(y), 1e-18)) if log_y else y
        return top + (ymax - value) / (ymax - ymin) * plot_h

    colors = ["#2563eb", "#dc2626", "#059669", "#7c3aed", "#d97706", "#0891b2"]
    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{title}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{left}" y="34" font-size="24" font-family="sans-serif" font-weight="700">{title}</text>',
    ]

    for i in range(6):
        y = top + i * plot_h / 5
        value = ymax - i * (ymax - ymin) / 5
        label = f"1e{value:.1f}" if log_y else f"{value:.5g}"
        chunks.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e5e7eb"/>')
        chunks.append(f'<text x="{left - 10}" y="{y + 5:.1f}" text-anchor="end" font-size="13" font-family="sans-serif" fill="#374151">{label}</text>')

    for i in range(6):
        x = left + i * plot_w / 5
        value = xmin + i * (xmax - xmin) / 5
        chunks.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" stroke="#f3f4f6"/>')
        chunks.append(f'<text x="{x:.1f}" y="{top + plot_h + 25}" text-anchor="middle" font-size="13" font-family="sans-serif" fill="#374151">{value:.4g}</text>')

    chunks.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#111827" stroke-width="1.5"/>')
    chunks.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#111827" stroke-width="1.5"/>')

    for idx, (label, values) in enumerate(series):
        points = []
        for x, y in zip(x_values, values):
            if math.isfinite(x) and math.isfinite(y):
                points.append(f"{sx(x):.2f},{sy(y):.2f}")
        if points:
            color = colors[idx % len(colors)]
            point_text = " ".join(points)
            chunks.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="{point_text}"/>')
            lx = left + 10 + (idx % 3) * 245
            ly = top + 20 + (idx // 3) * 24
            chunks.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 30}" y2="{ly}" stroke="{color}" stroke-width="3"/>')
            chunks.append(f'<text x="{lx + 38}" y="{ly + 5}" font-size="13" font-family="sans-serif" fill="#111827">{label}</text>')

    chunks.append(f'<text x="{left + plot_w / 2:.1f}" y="{height - 22}" text-anchor="middle" font-size="15" font-family="sans-serif">{x_label}</text>')
    chunks.append(f'<text x="22" y="{top + plot_h / 2:.1f}" text-anchor="middle" font-size="15" font-family="sans-serif" transform="rotate(-90 22 {top + plot_h / 2:.1f})">{y_label}</text>')
    chunks.append('</svg>')

    path.write_text("\n".join(chunks), encoding="utf-8")


def generate(taylor_csv, similarity_csv):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    tg = read_csv(taylor_csv)
    sim = read_csv(similarity_csv)

    times = [f(r, "time") for r in tg]
    svg_line_plot(
        FIG_DIR / "taylor_green_energy.svg",
        "Taylor–Green 基準計算: エネルギーとエンストロフィー",
        times,
        [
            ("運動エネルギー", [f(r, "kinetic_energy") for r in tg]),
            ("エンストロフィー", [f(r, "enstrophy") for r in tg]),
            ("粘性散逸", [f(r, "viscous_dissipation") for r in tg]),
        ],
        "時間 t",
        "診断量",
    )

    svg_line_plot(
        FIG_DIR / "taylor_green_errors.svg",
        "Taylor–Green 基準計算: 数値誤差診断",
        times,
        [
            ("発散 L2", [f(r, "divergence_l2") for r in tg]),
            ("PDE 離散残差 L2", [f(r, "pde_residual_l2") for r in tg]),
            ("射影後非線形項 L2", [f(r, "projected_nonlinear_l2") for r in tg]),
        ],
        "時間 t",
        "絶対値（対数表示）",
        log_y=True,
    )

    tau_axis = [-math.log10(f(r, "tau")) for r in sim]
    svg_line_plot(
        FIG_DIR / "similarity_coordinate_error.svg",
        "類似座標ソルバ: 論文恒等式の再構成誤差",
        tau_axis,
        [
            ("q 相対誤差", [f(r, "relative_q_error") for r in sim]),
            ("eta 絶対誤差", [f(r, "eta_error") for r in sim]),
        ],
        "-log10(tau)",
        "誤差（対数表示）",
        log_y=True,
    )

    last = tg[-1]
    max_div = max(f(r, "divergence_l2") for r in tg if math.isfinite(f(r, "divergence_l2")))
    max_cfl = max(f(r, "cfl") for r in tg if math.isfinite(f(r, "cfl")))
    max_q_err = max(f(r, "relative_q_error") for r in sim if math.isfinite(f(r, "relative_q_error")))
    max_eta_err = max(f(r, "eta_error") for r in sim if math.isfinite(f(r, "eta_error")))
    sha = os.environ.get("GITHUB_SHA", "local")[:12]

    report = f"""# 最新の自動計算レポート

このレポートは GitHub Actions で自動生成されます。対象コミット: `{sha}`

## 実施した計算

- FFTW 版 Taylor–Green 基準計算
- OpenAI 論文の式 (3.2)/(4.1) に対応する類似座標 `q, eta, X` の再構成テスト
- Lemma 4.1, 式 (4.2) の `T_b`, `Z_b` と物理座標有限差分の照合
- 式 (4.6)–(4.7) の `A_X`, `V0`, `Pi` 数値部品の manufactured-profile 回帰
- Fourier 微分・非圧縮射影・2/3 dealiasing・粘性減衰・Taylor–Green 回帰

## Taylor–Green 基準計算の要約

- 最終時刻: `{f(last, 'time'):.6g}`
- 最終運動エネルギー: `{f(last, 'kinetic_energy'):.12g}`
- 最終エンストロフィー: `{f(last, 'enstrophy'):.12g}`
- 最大発散誤差 L2: `{max_div:.3e}`
- 最大 CFL: `{max_cfl:.6g}`

![Taylor–Green energy](../results/figures/taylor_green_energy.svg)

![Taylor–Green errors](../results/figures/taylor_green_errors.svg)

## 類似座標ソルバの要約

論文の恒等式

`tau = q (1 - eta^2)`

`z = q^(1/2-h) eta`

を使って既知の `(q, eta)` から `(tau, z)` を作り、数値ソルバで `q` と `eta` を逆算しています。

- 最大 `q` 相対誤差: `{max_q_err:.3e}`
- 最大 `eta` 絶対誤差: `{max_eta_err:.3e}`

![Similarity coordinate errors](../results/figures/similarity_coordinate_error.svg)

## 解釈

論文固有の座標変換、Lemma 4.1 の微分作用素、式 (4.6)–(4.7) の radial average・非圧縮 radial flux・pressure integration の汎用数値部品まで実装しています。式 (4.8)–(4.11) の leading stress は `docs/leading-stress.md` に抽出済みです。

次の段階は、leading stress の数値部品を manufactured profile で検証し、その後 Theorem 4.6 と Appendices A/B の constructive profile を数値化することです。exact profile `E,U` を任意の surrogate で置き換えて OpenAI 構成と呼ぶことはしません。

CSV の生データは `results/reference/` に保存されています。
"""
    (REPORT_DIR / "latest.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: generate_artifacts.py <taylor_green.csv> <similarity_coordinates.csv>", file=sys.stderr)
        sys.exit(2)
    generate(Path(sys.argv[1]), Path(sys.argv[2]))
