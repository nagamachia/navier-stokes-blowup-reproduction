#!/usr/bin/env python3

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures"
REPORT_DIR = ROOT / "reports"


def font(size):
    for name in ("DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            pass
    return ImageFont.load_default()


def color(value, vmax):
    x = max(0.0, min(1.0, value / max(vmax, 1e-15)))
    # dark blue -> cyan -> yellow; compact and GIF-friendly
    if x < 0.5:
        a = 2.0 * x
        return (int(20 + 20 * a), int(35 + 150 * a), int(110 + 110 * a))
    a = 2.0 * (x - 0.5)
    return (int(40 + 215 * a), int(185 + 65 * a), int(220 - 180 * a))


def load(path):
    by_time_plane = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            t = float(row["time"])
            p = int(row["plane"])
            by_time_plane[(t, p)].append({k: float(v) if k not in ("plane",) else int(v) for k, v in row.items()})
    times = sorted({k[0] for k in by_time_plane})
    if not times:
        raise RuntimeError("snapshot CSV is empty")
    return by_time_plane, times


def axis_values(rows, plane):
    if plane == 0:
        return sorted({r["x"] for r in rows}), sorted({r["y"] for r in rows})
    if plane == 1:
        return sorted({r["x"] for r in rows}), sorted({r["z"] for r in rows})
    return sorted({r["y"] for r in rows}), sorted({r["z"] for r in rows})


def coords(row, plane):
    if plane == 0:
        return row["x"], row["y"]
    if plane == 1:
        return row["x"], row["z"]
    return row["y"], row["z"]


def frame_image(data, t, vmax):
    width, height = 1040, 390
    panel = 285
    top = 76
    lefts = [40, 377, 714]
    labels = ["xy slice (z=pi)", "xz slice (y=pi)", "yz slice (x=pi)"]
    im = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(im)
    d.text((40, 20), "3D Taylor-Green vortex: vorticity magnitude", fill=(20, 20, 20), font=font(24))
    d.text((760, 25), f"t = {t:.3f}", fill=(20, 20, 20), font=font(18))

    for plane in range(3):
        rows = data[(t, plane)]
        xs, ys = axis_values(rows, plane)
        xi = {v: i for i, v in enumerate(xs)}
        yi = {v: i for i, v in enumerate(ys)}
        cw = panel / max(len(xs), 1)
        ch = panel / max(len(ys), 1)
        ox = lefts[plane]
        for r in rows:
            x, y = coords(r, plane)
            i = xi[x]
            j = yi[y]
            x0 = ox + i * cw
            y0 = top + (len(ys) - 1 - j) * ch
            d.rectangle((x0, y0, x0 + cw + 1, y0 + ch + 1), fill=color(r["omega_mag"], vmax))
        d.rectangle((ox, top, ox + panel, top + panel), outline=(20, 20, 20), width=2)
        d.text((ox, top - 25), labels[plane], fill=(20, 20, 20), font=font(16))

    d.text((40, 365), f"shared scale: |omega| = 0 ... {vmax:.4g}", fill=(70, 70, 70), font=font(14))
    return im


def write_summary_svg(path, data, times):
    maxima = []
    for t in times:
        vals = [r["omega_mag"] for p in range(3) for r in data[(t, p)]]
        maxima.append(max(vals))
    w, h = 860, 420
    l, r, top, b = 85, 30, 55, 65
    pw, ph = w-l-r, h-top-b
    xmin, xmax = times[0], times[-1]
    ymin, ymax = 0.0, max(maxima) * 1.08
    sx = lambda x: l + (x-xmin)/max(xmax-xmin,1e-15)*pw
    sy = lambda y: top + (ymax-y)/max(ymax-ymin,1e-15)*ph
    pts = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x,y in zip(times,maxima))
    c = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}">', '<rect width="100%" height="100%" fill="white"/>',
         f'<text x="{l}" y="32" font-family="sans-serif" font-size="22" font-weight="700">3D Taylor–Green: max vorticity on central slices</text>']
    for q in range(6):
        y=top+q*ph/5; val=ymax-q*ymax/5
        c.append(f'<line x1="{l}" y1="{y:.2f}" x2="{l+pw}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        c.append(f'<text x="{l-10}" y="{y+5:.2f}" text-anchor="end" font-family="sans-serif" font-size="12">{val:.3g}</text>')
    c.append(f'<polyline points="{pts}" fill="none" stroke="#2563eb" stroke-width="2.5"/>')
    c.append(f'<line x1="{l}" y1="{top+ph}" x2="{l+pw}" y2="{top+ph}" stroke="#111827"/>')
    c.append(f'<line x1="{l}" y1="{top}" x2="{l}" y2="{top+ph}" stroke="#111827"/>')
    c.append(f'<text x="{l+pw/2}" y="{h-18}" text-anchor="middle" font-family="sans-serif" font-size="14">time t</text>')
    c.append('</svg>')
    path.write_text("\n".join(c), encoding="utf-8")
    return maxima


def main():
    if len(sys.argv) != 2:
        print("usage: generate_taylor_green_3d_gif.py <3d_snapshot.csv>", file=sys.stderr)
        raise SystemExit(2)
    source = Path(sys.argv[1])
    data, times = load(source)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    vmax = max(r["omega_mag"] for rows in data.values() for r in rows)
    frames = [frame_image(data, t, vmax) for t in times]
    gif_path = FIG_DIR / "taylor_green_3d_evolution.gif"
    frames[0].save(gif_path, save_all=True, append_images=frames[1:], duration=140, loop=0, optimize=False)
    maxima = write_summary_svg(FIG_DIR / "taylor_green_3d_max_vorticity.svg", data, times)
    report = f'''# 3D Taylor–Green 時間発展

標準的な3次元 Taylor–Green 初期条件を FFTW 擬スペクトル法 + RK4 で時間発展させ、中央 `xy / xz / yz` 断面の渦度強度を可視化しています。

## GIF アニメーション

![3D Taylor–Green evolution](../results/figures/taylor_green_3d_evolution.gif)

全 {len(times)} フレーム、`t={times[0]:.3f}` から `t={times[-1]:.3f}` までを共通色スケールで表示します。

## 最大渦度の推移

![max vorticity](../results/figures/taylor_green_3d_max_vorticity.svg)

- 初期中央断面最大渦度: `{maxima[0]:.8g}`
- 最終中央断面最大渦度: `{maxima[-1]:.8g}`
- 生データ: `results/reference/taylor_green_3d_snapshots_n24.csv`

この3D Taylor–Greenはソルバ・可視化基盤の検証用です。OpenAI論文の blow-up profile そのものではありません。
'''
    (REPORT_DIR / "taylor-green-3d.md").write_text(report, encoding="utf-8")
    print(f"wrote {gif_path} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
