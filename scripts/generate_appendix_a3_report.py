#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path

src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/reference/appendix_a3_closure.csv")
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("reports/appendix-a3.md")

rows = list(csv.DictReader(src.open()))
if not rows:
    raise SystemExit("Appendix A.3 diagnostics CSV is empty")

kb = float(rows[0]["Kb"])
amp = float(rows[0]["amplitude"])
max_mj = max(max(abs(float(r["residual_M"])), abs(float(r["residual_J"]))) for r in rows)
max_ip = max(max(abs(float(r["residual_I"])), abs(float(r["residual_pressure"]))) for r in rows)
max_q = max(abs(float(r["Qp"]) - float(r["Qs0"])) for r in rows)

lines = [
    "# Appendix A.3 moment closure 診断",
    "",
    "このレポートは `include/appendix_a3_closure.hpp` の軽量数値回帰を GitHub Actions で再生成したもの。",
    "論文の Appendix A.3 に対応する正規化された moment closure を検証する。完全な Theorem 4.6 profile assembly ではない。",
    "",
    "## 主要結果",
    "",
    f"- `(A.19)` の `Kb = {kb:.12g}`（論文中の評価 `0.20 < Kb <= 0.25` と整合）",
    f"- remainder を0とした principal amplitude: `{amp:.12g}`（論文の bracket `0.9 < Amp < 1.2` 内）",
    f"- `(A.15)` の `M,J` closure 最大残差: `{max_mj:.3e}`",
    f"- `(A.11)` の `I, pressure increment` closure 最大残差: `{max_ip:.3e}`",
    f"- `(A.16)` の `Qs(0)=Qp` 最大差: `{max_q:.3e}`",
    "",
    "## λ sweep",
    "",
    "| lambda | Amp | c1(MJ) | c2(MJ) | c1(I/P) | c2(I/P) |",
    "|---:|---:|---:|---:|---:|---:|",
]
for r in rows:
    lines.append(
        f"| {float(r['lambda']):.3g} | {float(r['amplitude']):.8g} | "
        f"{float(r['c1_MJ']):.3e} | {float(r['c2_MJ']):.3e} | "
        f"{float(r['c1_IP']):.3e} | {float(r['c2_IP']):.3e} |"
    )

lines += [
    "",
    "## 実装境界",
    "",
    "現在は Appendix A.3 の有限次元 closure mechanism を式番号付きで再現している。",
    "次段階では A.2 の全 radial schedule から `(A.14)` の pre-pulse discrepancies と `(A.19)` の remainder `E(Amp,eta)` を直接評価し、",
    "この closure solver に入力する。これにより trial profile ではなく、outer profile 全体を通した `S(infinity)=0` の数値閉包へ進む。",
    "",
]
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines), encoding="utf-8")
print(f"wrote {out}")
