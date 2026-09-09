#!/usr/bin/env python3

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures"
REF_DIR = ROOT / "results" / "reference"
REPORT_DIR = ROOT / "reports"
TWO_PI = 2.0 * math.pi


def color_diverging(value, vmax):
    vmax = max(vmax, 1e-30)
    x = max(-1.0, min(1.0, value / vmax))
    if x >= 0:
        r = 245
        g = int(245 - 150 * x)
        b = int(245 - 170 * x)
    else:
        x = -x
        r = int(245 - 170 * x)
        g = int(245 - 135 * x)
        b = 245
    return f"rgb({r},{g},{b})"


def color_error(value, vmax):
    vmax = max(vmax, 1e-30)
    x = max(0.0, min(1.0, value / vmax))
    r = int(245 - 35 * x)
    g = int(245 - 180 * x)
    b = int(245 - 205 * x)
    return f"rgb({r},{g},{b})"


def read_snapshots(path):
    groups = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            t = float(row["time"])
            groups[t].append({
                "i": int(row["i"]),
                "j": int(row["j"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "u": float(row["u"]),
                "v": float(row["v"]),
                "omega_z": float(row["omega_z"]),
                "exact_u": float(row["exact_u"]),
                "exact_v": float(row["exact_v"]),
                "velocity_error": float(row["velocity_error"]),
            })
    if not groups:
        raise RuntimeError(f"no snapshot rows in {path}")
    return dict(sorted(groups.items()))


def grid_n(rows):
    return 1 + max(max(r["i"], r["j"]) for r in rows)


def field_panel_svg(path, snapshots, key, title, subtitle, error_mode=False, arrows=False):
    times = list(snapshots.keys())
    panel = 250
    margin = 48
    gap = 16
    width = margin * 2 + len(times) * panel + (len(times) - 1) * gap
    height = 380
    values = [abs(r[key]) for rows in snapshots.values() for r in rows]
    vmax = max(values) if values else 1.0
    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="48" y="34" font-family="sans-serif" font-size="24" font-weight="700">{title}</text>',
        f'<text x="48" y="60" font-family="sans-serif" font-size="14">{subtitle}</text>',
    ]

    for p, t in enumerate(times):
        rows = snapshots[t]
        n = grid_n(rows)
        ox = margin + p * (panel + gap)
        oy = 84
        cell = panel / n
        rowmap = {(r["i"], r["j"]): r for r in rows}
        for j in range(n):
            for i in range(n):
                r = rowmap[(i, j)]
                value = r[key]
                fill = color_error(value, vmax) if error_mode else color_diverging(value, vmax)
                chunks.append(
                    f'<rect x="{ox + i*cell:.2f}" y="{oy + (n-1-j)*cell:.2f}" '
                    f'width="{cell+.25:.2f}" height="{cell+.25:.2f}" fill="{fill}"/>'
                )
        if arrows:
            stride = max(1, n // 8)
            scale = 15.0
            chunks.append('<g stroke="#111827" stroke-width="1.1" fill="none" opacity="0.82">')
            for j in range(0, n, stride):
                for i in range(0, n, stride):
                    r = rowmap[(i, j)]
                    u, v = r["u"], r["v"]
                    cx = ox + panel * (i + 0.5) / n
                    cy = oy + panel * (1.0 - (j + 0.5) / n)
                    dx, dy = scale * u, -scale * v
                    chunks.append(f'<line x1="{cx:.2f}" y1="{cy:.2f}" x2="{cx+dx:.2f}" y2="{cy+dy:.2f}"/>')
            chunks.append('</g>')
        chunks.append(f'<rect x="{ox}" y="{oy}" width="{panel}" height="{panel}" fill="none" stroke="#111827"/>')
        chunks.append(f'<text x="{ox+panel/2:.1f}" y="{oy+panel+23}" text-anchor="middle" font-family="sans-serif" font-size="14">t = {t:.5f}</text>')
    chunks.append('</svg>')
    path.write_text('\n'.join(chunks), encoding='utf-8')


def error_history_svg(path, snapshots):
    times = list(snapshots.keys())
    max_err = [max(r["velocity_error"] for r in snapshots[t]) for t in times]
    rms_err = [math.sqrt(sum(r["velocity_error"]**2 for r in snapshots[t]) / len(snapshots[t])) for t in times]
    width, height = 860, 480
    left, right, top, bottom = 90, 30, 70, 70
    pw, ph = width-left-right, height-top-bottom
    ymax = max(max_err) * 1.15 if max(max_err) > 0 else 1.0
    xmin, xmax = min(times), max(times)
    def sx(x): return left + (x-xmin)/(xmax-xmin if xmax>xmin else 1.0)*pw
    def sy(y): return top + (ymax-y)/ymax*ph
    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="60" y="34" font-family="sans-serif" font-size="24" font-weight="700">Taylor–Green 数値解と解析解の速度誤差</text>',
        f'<line x1="{left}" y1="{top+ph}" x2="{left+pw}" y2="{top+ph}" stroke="#111827"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+ph}" stroke="#111827"/>',
    ]
    for frac in range(6):
        y = ymax * frac / 5
        py = sy(y)
        chunks.append(f'<line x1="{left}" y1="{py:.2f}" x2="{left+pw}" y2="{py:.2f}" stroke="#e5e7eb"/>')
        chunks.append(f'<text x="{left-8}" y="{py+5:.2f}" text-anchor="end" font-family="sans-serif" font-size="12">{y:.2e}</text>')
    for label, vals, dash in [("最大誤差", max_err, ""), ("RMS誤差", rms_err, "6 4")]:
        pts = ' '.join(f'{sx(t):.2f},{sy(v):.2f}' for t, v in zip(times, vals))
        chunks.append(f'<polyline points="{pts}" fill="none" stroke="#111827" stroke-width="2.4" stroke-dasharray="{dash}"/>')
        chunks.append(f'<text x="{left+pw-130}" y="{top+22+(22 if dash else 0)}" font-family="sans-serif" font-size="13">{label}</text>')
    chunks.append(f'<text x="{left+pw/2}" y="{height-20}" text-anchor="middle" font-family="sans-serif" font-size="14">時間 t</text>')
    chunks.append('</svg>')
    path.write_text('\n'.join(chunks), encoding='utf-8')


def write_numerical_vtk(path, rows, nz=8):
    n = grid_n(rows)
    rowmap = {(r["i"], r["j"]): r for r in rows}
    spacing = TWO_PI / n
    with path.open('w', encoding='utf-8') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Numerical Taylor-Green snapshot extruded in z\n')
        f.write('ASCII\nDATASET STRUCTURED_POINTS\n')
        f.write(f'DIMENSIONS {n} {n} {nz}\nORIGIN 0 0 0\n')
        f.write(f'SPACING {spacing:.17g} {spacing:.17g} {spacing:.17g}\n')
        f.write(f'POINT_DATA {n*n*nz}\n')
        f.write('VECTORS velocity double\n')
        for _k in range(nz):
            for j in range(n):
                for i in range(n):
                    r = rowmap[(i, j)]
                    f.write(f'{r["u"]:.17g} {r["v"]:.17g} 0\n')
        f.write('SCALARS omega_z double 1\nLOOKUP_TABLE default\n')
        for _k in range(nz):
            for j in range(n):
                for i in range(n):
                    f.write(f'{rowmap[(i,j)]["omega_z"]:.17g}\n')
        f.write('SCALARS velocity_error double 1\nLOOKUP_TABLE default\n')
        for _k in range(nz):
            for j in range(n):
                for i in range(n):
                    f.write(f'{rowmap[(i,j)]["velocity_error"]:.17g}\n')


def write_report(path, snapshots, csv_name, vtk_name):
    times = list(snapshots.keys())
    final_rows = snapshots[times[-1]]
    max_err = max(r["velocity_error"] for r in final_rows)
    rms_err = math.sqrt(sum(r["velocity_error"]**2 for r in final_rows) / len(final_rows))
    path.write_text(f'''# Taylor–Green 数値流れ可視化レポート

このレポートは解析式だけではなく、FFTW を使った数値時間積分が出力した格子スナップショットから生成しています。

## 数値渦度の時間発展

![numerical evolution](../results/figures/taylor_green_numerical_evolution.svg)

## 数値解と解析解の速度誤差分布

![numerical error field](../results/figures/taylor_green_velocity_error_field.svg)

## 誤差の時間推移

![numerical error history](../results/figures/taylor_green_velocity_error_history.svg)

最終スナップショットの最大速度誤差は `{max_err:.3e}`、RMS 速度誤差は `{rms_err:.3e}` です。

## 3D 可視化

`results/reference/{vtk_name}` は最終数値格子場を z 方向へ埋め込んだ ParaView 用 VTK です。`velocity`, `omega_z`, `velocity_error` を含みます。

ParaView では Slice / Glyph / Stream Tracer を使うと流れ方向を確認できます。CFD データに対して ParaView は流線や方向付き Glyph を使った可視化をサポートしています。

生データは `results/reference/{csv_name}` に保存しています。
''', encoding='utf-8')


def main():
    if len(sys.argv) != 2:
        print('usage: generate_taylor_green_evolution.py <numerical_snapshots.csv>', file=sys.stderr)
        raise SystemExit(2)
    csv_path = Path(sys.argv[1])
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    snapshots = read_snapshots(csv_path)
    field_panel_svg(
        FIG_DIR / 'taylor_green_numerical_evolution.svg', snapshots, 'omega_z',
        'Taylor–Green 数値渦度の時間発展',
        'FFTW + RK4 の格子場。背景: ω_z、矢印: 数値速度。', arrows=True)
    field_panel_svg(
        FIG_DIR / 'taylor_green_velocity_error_field.svg', snapshots, 'velocity_error',
        'Taylor–Green 数値解と解析解の速度誤差分布',
        '各格子点の |u_num - u_exact|。全時刻で共通スケール。', error_mode=True)
    error_history_svg(FIG_DIR / 'taylor_green_velocity_error_history.svg', snapshots)
    final_rows = snapshots[list(snapshots.keys())[-1]]
    vtk_name = 'taylor_green_numerical_n32.vtk'
    write_numerical_vtk(REF_DIR / vtk_name, final_rows)
    write_report(REPORT_DIR / 'flow-visualization.md', snapshots, csv_path.name, vtk_name)


if __name__ == '__main__':
    main()
