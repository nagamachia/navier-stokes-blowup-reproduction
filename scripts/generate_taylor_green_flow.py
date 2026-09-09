#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures"
REPORT = ROOT / "reports" / "latest.md"
TWO_PI = 2.0 * math.pi


def read_last_row(csv_path):
    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"CSV is empty: {csv_path}")
    return rows[-1]


def velocity(x, y, time, viscosity):
    decay = math.exp(-2.0 * viscosity * time)
    u = decay * math.sin(x) * math.cos(y)
    v = -decay * math.cos(x) * math.sin(y)
    return u, v


def vorticity(x, y, time, viscosity):
    decay = math.exp(-2.0 * viscosity * time)
    return 2.0 * decay * math.sin(x) * math.sin(y)


def color_speed(value, vmax):
    a = 0.0 if vmax <= 0 else max(0.0, min(1.0, value / vmax))
    r = int(245 - 80 * a)
    g = int(248 - 105 * a)
    b = int(255 - 10 * a)
    return f"rgb({r},{g},{b})"


def color_vorticity(value, vmax):
    a = 0.0 if vmax <= 0 else max(-1.0, min(1.0, value / vmax))
    if a >= 0:
        r = int(245 - 15 * a)
        g = int(245 - 145 * a)
        b = int(245 - 175 * a)
    else:
        aa = -a
        r = int(245 - 175 * aa)
        g = int(245 - 100 * aa)
        b = int(245 - 10 * aa)
    return f"rgb({r},{g},{b})"


def base_svg(title, subtitle=""):
    width, height = 760, 720
    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{title}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="70" y="34" font-size="23" font-family="sans-serif" font-weight="700">{title}</text>',
    ]
    if subtitle:
        chunks.append(f'<text x="70" y="56" font-size="13" font-family="sans-serif" fill="#4b5563">{subtitle}</text>')
    return chunks, width, height


def plot_geometry():
    left, top, size = 82.0, 82.0, 590.0
    return left, top, size


def sx(x):
    left, _, size = plot_geometry()
    return left + (x / TWO_PI) * size


def sy(y):
    _, top, size = plot_geometry()
    return top + (1.0 - y / TWO_PI) * size


def add_axes(chunks):
    left, top, size = plot_geometry()
    chunks.append(f'<rect x="{left}" y="{top}" width="{size}" height="{size}" fill="none" stroke="#111827" stroke-width="1.4"/>')
    for value, label in [(0.0, "0"), (math.pi / 2, "π/2"), (math.pi, "π"), (3 * math.pi / 2, "3π/2"), (TWO_PI, "2π")]:
        x = sx(value)
        y = sy(value)
        chunks.append(f'<line x1="{x:.2f}" y1="{top + size}" x2="{x:.2f}" y2="{top + size + 6}" stroke="#111827"/>')
        chunks.append(f'<text x="{x:.2f}" y="{top + size + 24}" text-anchor="middle" font-size="12" font-family="sans-serif">{label}</text>')
        chunks.append(f'<line x1="{left - 6}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#111827"/>')
        chunks.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-size="12" font-family="sans-serif">{label}</text>')
    chunks.append(f'<text x="{left + size / 2}" y="{top + size + 45}" text-anchor="middle" font-size="14" font-family="sans-serif">x</text>')
    chunks.append(f'<text x="28" y="{top + size / 2}" text-anchor="middle" font-size="14" font-family="sans-serif" transform="rotate(-90 28 {top + size / 2})">y</text>')


def write_velocity_field(path, time, viscosity):
    chunks, _, _ = base_svg(
        "Taylor–Green 速度場",
        f"z=0 断面, t={time:.3g}, ν={viscosity:.3g}。背景は速度の大きさ、矢印は流向。",
    )
    left, top, size = plot_geometry()
    cells = 40
    vmax = math.exp(-2.0 * viscosity * time)
    cell = size / cells
    for i in range(cells):
        for j in range(cells):
            x = (i + 0.5) * TWO_PI / cells
            y = (j + 0.5) * TWO_PI / cells
            u, v = velocity(x, y, time, viscosity)
            speed = math.hypot(u, v)
            chunks.append(
                f'<rect x="{left + i * cell:.2f}" y="{top + (cells - 1 - j) * cell:.2f}" width="{cell + 0.5:.2f}" height="{cell + 0.5:.2f}" fill="{color_speed(speed, vmax)}" stroke="none"/>'
            )

    chunks.append('<defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><polygon points="0 0, 7 3.5, 0 7" fill="#111827"/></marker></defs>')
    n = 17
    arrow_scale = 18.0
    for i in range(n):
        for j in range(n):
            x = (i + 0.5) * TWO_PI / n
            y = (j + 0.5) * TWO_PI / n
            u, v = velocity(x, y, time, viscosity)
            speed = math.hypot(u, v)
            if speed < 1e-10:
                continue
            length = arrow_scale * (0.35 + 0.65 * speed / vmax)
            dx = length * u / speed
            dy = -length * v / speed
            px, py = sx(x), sy(y)
            chunks.append(f'<line x1="{px - 0.5 * dx:.2f}" y1="{py - 0.5 * dy:.2f}" x2="{px + 0.5 * dx:.2f}" y2="{py + 0.5 * dy:.2f}" stroke="#111827" stroke-width="1.2" marker-end="url(#arrow)"/>')

    add_axes(chunks)
    chunks.append('</svg>')
    path.write_text("\n".join(chunks), encoding="utf-8")


def rk4_streamline(x, y, time, viscosity, direction, steps=520, ds=0.035):
    points = []

    def unit_field(px, py):
        u, v = velocity(px, py, time, viscosity)
        mag = math.hypot(u, v)
        if mag < 1e-8:
            return 0.0, 0.0
        return direction * u / mag, direction * v / mag

    for _ in range(steps):
        points.append((x, y))
        k1x, k1y = unit_field(x, y)
        if k1x == 0.0 and k1y == 0.0:
            break
        k2x, k2y = unit_field((x + 0.5 * ds * k1x) % TWO_PI, (y + 0.5 * ds * k1y) % TWO_PI)
        k3x, k3y = unit_field((x + 0.5 * ds * k2x) % TWO_PI, (y + 0.5 * ds * k2y) % TWO_PI)
        k4x, k4y = unit_field((x + ds * k3x) % TWO_PI, (y + ds * k3y) % TWO_PI)
        nx = (x + ds * (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0) % TWO_PI
        ny = (y + ds * (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0) % TWO_PI
        if math.hypot(nx - x, ny - y) > 1.0:
            break
        x, y = nx, ny
    return points


def write_streamlines(path, time, viscosity):
    chunks, _, _ = base_svg(
        "Taylor–Green 流線",
        "閉じたセル状の循環が交互に並ぶ。流線は速度場を数値積分して描画。",
    )
    left, top, size = plot_geometry()
    chunks.append(f'<rect x="{left}" y="{top}" width="{size}" height="{size}" fill="#f8fafc"/>')
    seeds = []
    for cx in [math.pi / 2, 3 * math.pi / 2]:
        for cy in [math.pi / 2, 3 * math.pi / 2]:
            for radius in [0.28, 0.52, 0.78, 1.02]:
                seeds.append((cx + radius, cy))

    for seed_x, seed_y in seeds:
        backward = rk4_streamline(seed_x % TWO_PI, seed_y % TWO_PI, time, viscosity, -1.0, steps=250)
        forward = rk4_streamline(seed_x % TWO_PI, seed_y % TWO_PI, time, viscosity, 1.0, steps=250)
        pts = list(reversed(backward)) + forward[1:]
        if len(pts) < 2:
            continue
        rendered = []
        previous = None
        for x, y in pts:
            p = (sx(x), sy(y))
            if previous is not None and math.hypot(p[0] - previous[0], p[1] - previous[1]) > 80:
                if len(rendered) > 1:
                    chunks.append(f'<polyline fill="none" stroke="#1f2937" stroke-width="1.25" opacity="0.82" points="{" ".join(rendered)}"/>')
                rendered = []
            rendered.append(f"{p[0]:.2f},{p[1]:.2f}")
            previous = p
        if len(rendered) > 1:
            chunks.append(f'<polyline fill="none" stroke="#1f2937" stroke-width="1.25" opacity="0.82" points="{" ".join(rendered)}"/>')

    for x in [math.pi, 2 * math.pi]:
        chunks.append(f'<line x1="{sx(x):.2f}" y1="{top}" x2="{sx(x):.2f}" y2="{top + size}" stroke="#94a3b8" stroke-dasharray="4 5"/>')
    for y in [math.pi, 2 * math.pi]:
        chunks.append(f'<line x1="{left}" y1="{sy(y):.2f}" x2="{left + size}" y2="{sy(y):.2f}" stroke="#94a3b8" stroke-dasharray="4 5"/>')

    add_axes(chunks)
    chunks.append('</svg>')
    path.write_text("\n".join(chunks), encoding="utf-8")


def write_vorticity(path, time, viscosity):
    chunks, _, _ = base_svg(
        "Taylor–Green 渦度 ωz",
        "暖色と寒色は逆向きの回転。ωz = 2 exp(-2νt) sin(x) sin(y)。",
    )
    left, top, size = plot_geometry()
    cells = 64
    vmax = 2.0 * math.exp(-2.0 * viscosity * time)
    cell = size / cells
    for i in range(cells):
        for j in range(cells):
            x = (i + 0.5) * TWO_PI / cells
            y = (j + 0.5) * TWO_PI / cells
            omega = vorticity(x, y, time, viscosity)
            chunks.append(
                f'<rect x="{left + i * cell:.2f}" y="{top + (cells - 1 - j) * cell:.2f}" width="{cell + 0.4:.2f}" height="{cell + 0.4:.2f}" fill="{color_vorticity(omega, vmax)}" stroke="none"/>'
            )
    add_axes(chunks)
    chunks.append('</svg>')
    path.write_text("\n".join(chunks), encoding="utf-8")


def append_report(time, viscosity):
    if not REPORT.exists():
        return
    text = REPORT.read_text(encoding="utf-8")
    marker = "## Taylor–Green の流れの可視化"
    if marker in text:
        text = text.split(marker, 1)[0].rstrip() + "\n"
    section = f"""

## Taylor–Green の流れの可視化

基準計算と同じ `t={time:.3g}`, `nu={viscosity:.3g}` の `z=0` 断面を可視化しています。現在の Taylor–Green 回帰は z に依存しない2次元渦を3次元周期箱へ埋め込んだものなので、この断面が全 z で同じ形になります。

### 速度場

背景色は速度の大きさ、矢印は速度ベクトルです。

![Taylor–Green velocity field](../results/figures/taylor_green_velocity_field.svg)

### 流線

閉じた循環セルが交互に並び、隣接セルでは回転方向が反転します。

![Taylor–Green streamlines](../results/figures/taylor_green_streamlines.svg)

### 渦度

`omega_z` の符号で回転方向が分かります。暖色・寒色の4領域が交互に現れます。

![Taylor–Green vorticity](../results/figures/taylor_green_vorticity.svg)
"""
    REPORT.write_text(text.rstrip() + section + "\n", encoding="utf-8")


def main():
    if len(sys.argv) not in (2, 3):
        print("usage: generate_taylor_green_flow.py <taylor_green.csv> [viscosity]", file=sys.stderr)
        return 2
    last = read_last_row(Path(sys.argv[1]))
    time = float(last["time"])
    viscosity = float(sys.argv[2]) if len(sys.argv) == 3 else 0.1
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    write_velocity_field(FIG_DIR / "taylor_green_velocity_field.svg", time, viscosity)
    write_streamlines(FIG_DIR / "taylor_green_streamlines.svg", time, viscosity)
    write_vorticity(FIG_DIR / "taylor_green_vorticity.svg", time, viscosity)
    append_report(time, viscosity)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
