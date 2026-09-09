#!/usr/bin/env python3

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures"
REF_DIR = ROOT / "results" / "reference"
REPORT_DIR = ROOT / "reports"

TWO_PI = 2.0 * math.pi


def velocity(x, y, t, nu):
    a = math.exp(-2.0 * nu * t)
    u = a * math.sin(x) * math.cos(y)
    v = -a * math.cos(x) * math.sin(y)
    return u, v, 0.0


def omega_z(x, y, t, nu):
    return 2.0 * math.exp(-2.0 * nu * t) * math.sin(x) * math.sin(y)


def color_for(value, vmax):
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


def panel_vorticity_svg(path, times, nu):
    panel = 260
    margin = 54
    gap = 18
    width = margin * 2 + len(times) * panel + (len(times) - 1) * gap
    height = 390
    n = 40
    vmax = 2.0

    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="54" y="34" font-family="sans-serif" font-size="24" font-weight="700">Taylor–Green 渦度の時間発展（共通色スケール）</text>',
        '<text x="54" y="60" font-family="sans-serif" font-size="14">背景: ω_z、矢印: 速度方向。形は保たれ、粘性で振幅が減衰する。</text>',
    ]

    for p, t in enumerate(times):
        ox = margin + p * (panel + gap)
        oy = 88
        cell = panel / n
        for j in range(n):
            y = TWO_PI * (j + 0.5) / n
            for i in range(n):
                x = TWO_PI * (i + 0.5) / n
                wz = omega_z(x, y, t, nu)
                chunks.append(
                    f'<rect x="{ox + i * cell:.2f}" y="{oy + (n - 1 - j) * cell:.2f}" '
                    f'width="{cell + 0.3:.2f}" height="{cell + 0.3:.2f}" fill="{color_for(wz, vmax)}"/>'
                )

        # velocity arrows on a coarser grid
        m = 8
        scale = 14.0
        chunks.append('<g stroke="#111827" stroke-width="1.2" fill="none" opacity="0.85">')
        for j in range(m):
            y = TWO_PI * (j + 0.5) / m
            for i in range(m):
                x = TWO_PI * (i + 0.5) / m
                u, v, _ = velocity(x, y, t, nu)
                norm = math.hypot(u, v)
                if norm < 1e-10:
                    continue
                cx = ox + panel * (i + 0.5) / m
                cy = oy + panel * (1.0 - (j + 0.5) / m)
                dx = scale * u
                dy = -scale * v
                chunks.append(f'<line x1="{cx:.2f}" y1="{cy:.2f}" x2="{cx + dx:.2f}" y2="{cy + dy:.2f}"/>')
                ang = math.atan2(dy, dx)
                ah = 4.0
                for sgn in (-1, 1):
                    aa = ang + math.pi + sgn * 0.45
                    x2 = cx + dx + ah * math.cos(aa)
                    y2 = cy + dy + ah * math.sin(aa)
                    chunks.append(f'<line x1="{cx + dx:.2f}" y1="{cy + dy:.2f}" x2="{x2:.2f}" y2="{y2:.2f}"/>')
        chunks.append('</g>')
        chunks.append(f'<rect x="{ox}" y="{oy}" width="{panel}" height="{panel}" fill="none" stroke="#111827"/>')
        chunks.append(f'<text x="{ox + panel/2:.1f}" y="{oy + panel + 24}" text-anchor="middle" font-family="sans-serif" font-size="15">t = {t:.3f}</text>')

    chunks.append('</svg>')
    path.write_text('\n'.join(chunks), encoding='utf-8')


def stacked_slices_svg(path, t, nu):
    width, height = 920, 620
    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="48" y="38" font-family="sans-serif" font-size="24" font-weight="700">Taylor–Green の 3D 構造（2D 渦の z 方向埋め込み）</text>',
        '<text x="48" y="66" font-family="sans-serif" font-size="14">現在の基準解は z に依存しないため、同じ xy 渦パターンが z 方向に積層される。</text>',
    ]
    n = 18
    size = 270
    slices = 6
    for s in range(slices - 1, -1, -1):
        ox = 150 + s * 55
        oy = 240 - s * 28
        cell = size / n
        for j in range(n):
            y = TWO_PI * (j + 0.5) / n
            for i in range(n):
                x = TWO_PI * (i + 0.5) / n
                wz = omega_z(x, y, t, nu)
                chunks.append(
                    f'<rect x="{ox + i*cell:.2f}" y="{oy + (n-1-j)*cell:.2f}" width="{cell+.2:.2f}" height="{cell+.2:.2f}" '
                    f'fill="{color_for(wz, 2.0)}" opacity="0.82"/>'
                )
        chunks.append(f'<rect x="{ox}" y="{oy}" width="{size}" height="{size}" fill="none" stroke="#111827" stroke-width="1.2"/>')
        chunks.append(f'<text x="{ox + size + 10}" y="{oy + 14}" font-family="sans-serif" font-size="13">z slice {s}</text>')
    chunks.append('<text x="150" y="555" font-family="sans-serif" font-size="14">実際の3D爆発候補では、この z 不変性は失われるため、将来は等値面・3D流線へ置き換える。</text>')
    chunks.append('</svg>')
    path.write_text('\n'.join(chunks), encoding='utf-8')


def write_vtk(path, n, t, nu):
    # Legacy ASCII VTK: easy to inspect and opens directly in ParaView.
    with path.open('w', encoding='utf-8') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write(f'Taylor-Green reference t={t} nu={nu}\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_POINTS\n')
        f.write(f'DIMENSIONS {n} {n} {n}\n')
        f.write('ORIGIN 0 0 0\n')
        spacing = TWO_PI / n
        f.write(f'SPACING {spacing:.17g} {spacing:.17g} {spacing:.17g}\n')
        f.write(f'POINT_DATA {n*n*n}\n')
        f.write('VECTORS velocity double\n')
        for k in range(n):
            z = TWO_PI * k / n
            _ = z
            for j in range(n):
                y = TWO_PI * j / n
                for i in range(n):
                    x = TWO_PI * i / n
                    u, v, w = velocity(x, y, t, nu)
                    f.write(f'{u:.17g} {v:.17g} {w:.17g}\n')
        f.write('SCALARS omega_z double 1\n')
        f.write('LOOKUP_TABLE default\n')
        for k in range(n):
            for j in range(n):
                y = TWO_PI * j / n
                for i in range(n):
                    x = TWO_PI * i / n
                    f.write(f'{omega_z(x, y, t, nu):.17g}\n')


def write_report(path):
    path.write_text('''# Taylor–Green 可視化レポート

## 時間発展

5時刻を同じ色スケールで並べ、渦度の減衰を比較します。

![Taylor–Green evolution](../results/figures/taylor_green_evolution.svg)

## 3D 構造

現在の基準 Taylor–Green は 2D 解を 3D 周期箱へ埋め込んだものなので、`z` 方向には変化しません。下図は同じ `xy` 渦パターンが積層されることを示します。

![Taylor–Green stacked slices](../results/figures/taylor_green_stacked_slices.svg)

## ParaView 用データ

`results/reference/taylor_green_reference_n16.vtk` を ParaView で開くと、速度ベクトル `velocity` と渦度 `omega_z` を3次元で確認できます。

この VTK は現在、解析的 Taylor–Green 基準解から生成しています。次段階では FFTW ソルバ本体の実格子 snapshot を同じ形式で出力し、解析解との差も可視化します。
''', encoding='utf-8')


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    nu = 0.1
    times = [0.0, 0.025, 0.05, 0.075, 0.1]
    panel_vorticity_svg(FIG_DIR / 'taylor_green_evolution.svg', times, nu)
    stacked_slices_svg(FIG_DIR / 'taylor_green_stacked_slices.svg', 0.1, nu)
    write_vtk(REF_DIR / 'taylor_green_reference_n16.vtk', 16, 0.1, nu)
    write_report(REPORT_DIR / 'flow-visualization.md')


if __name__ == '__main__':
    main()
