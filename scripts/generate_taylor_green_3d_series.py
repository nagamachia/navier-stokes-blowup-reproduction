#!/usr/bin/env python3

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

TWO_PI = 2.0 * math.pi


def load_snapshots(path: Path):
    groups = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            groups[float(row["time"])].append(row)
    if not groups:
        raise RuntimeError(f"no rows in {path}")
    return dict(sorted(groups.items()))


def write_vtk(path: Path, rows):
    n = 1 + max(max(int(r["i"]), int(r["j"]), int(r["k"])) for r in rows)
    if len(rows) != n ** 3:
        raise RuntimeError(f"expected {n**3} rows, got {len(rows)}")
    spacing = TWO_PI / n
    rowmap = {(int(r["i"]), int(r["j"]), int(r["k"])): r for r in rows}
    with path.open("w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("3D Taylor-Green time-series snapshot\n")
        f.write("ASCII\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {n} {n} {n}\n")
        f.write("ORIGIN 0 0 0\n")
        f.write(f"SPACING {spacing:.17g} {spacing:.17g} {spacing:.17g}\n")
        f.write(f"POINT_DATA {n**3}\n")
        f.write("VECTORS velocity double\n")
        for k in range(n):
            for j in range(n):
                for i in range(n):
                    r = rowmap[(i, j, k)]
                    f.write(f'{r["u"]} {r["v"]} {r["w"]}\n')
        f.write("SCALARS omega_mag double 1\nLOOKUP_TABLE default\n")
        for k in range(n):
            for j in range(n):
                for i in range(n):
                    f.write(f'{rowmap[(i,j,k)]["omega_mag"]}\n')
    return n


def main():
    if len(sys.argv) != 3:
        print("usage: generate_taylor_green_3d_series.py <full_snapshot.csv> <output_dir>", file=sys.stderr)
        raise SystemExit(2)
    source = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)
    snapshots = load_snapshots(source)
    frames = []
    grid_n = None
    for index, (time, rows) in enumerate(snapshots.items()):
        name = f"frame_{index:03d}.vtk"
        n = write_vtk(output_dir / name, rows)
        if grid_n is None:
            grid_n = n
        elif grid_n != n:
            raise RuntimeError("grid size changed between frames")
        frames.append({"index": index, "time": time, "file": name})
    manifest = {
        "format": "legacy-vtk-structured-points-series-v1",
        "grid": [grid_n, grid_n, grid_n],
        "frame_count": len(frames),
        "frames": frames,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(frames)} VTK frames to {output_dir}")


if __name__ == "__main__":
    main()
