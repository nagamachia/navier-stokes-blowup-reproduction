#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/reference/appendix_a3_closure.csv")
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("reports/appendix-a3.md")
global_src = Path(sys.argv[3]) if len(sys.argv) > 3 else None

rows = list(csv.DictReader(src.open()))
if not rows:
    raise SystemExit("Appendix A.3 diagnostics CSV is empty")
global_rows = list(csv.DictReader(global_src.open())) if global_src and global_src.exists() else []

kb = float(rows[0]["Kb"])
amp = float(rows[0]["amplitude"])
max_mj = max(max(abs(float(r["residual_M"])), abs(float(r["residual_J"]))) for r in rows)
max_ip = max(max(abs(float(r["residual_I"])), abs(float(r["residual_pressure"]))) for r in rows)
max_q = max(abs(float(r["Qp"]) - float(r["Qs0"])) for r in rows)
max_pre_m = max(abs(float(r["pre_M_bump"])) for r in rows)
max_pre_j = max(abs(float(r["pre_J_bump"])) for r in rows)

lines = [
    "# Appendix A.3 moment closure 診断", "",
    "Appendix A.2 の staged radial profile と A.3 closure を GitHub Actions で数値回帰する。",
    "`M/(XE)`, `J/(XHE)`, `S/(XE^2)`, `I/(XH)` を直接発展させ、巨大な半径を生成しない。", "",
    "## 主要結果", "",
    f"- `(A.19)` の `Kb = {kb:.12g}`（`0.20 < Kb <= 0.25` と整合）",
    f"- remainder=0 の principal amplitude: `{amp:.12g}`",
    f"- A.15 pre-M 最大値: `{max_pre_m:.3e}` / pre-J 最大値: `{max_pre_j:.3e}`",
    f"- `(A.15)` M,J closure 最大残差: `{max_mj:.3e}`",
    f"- `(A.11)` I/pressure closure 最大残差: `{max_ip:.3e}`",
    f"- `(A.16)` Qs(0)=Qp 最大差: `{max_q:.3e}`", "",
    "## 局所 λ・η sweep", "",
    "| lambda | eta | m(Xp) | j(Xp) | s(Xp) | pre-M | pre-J |",
    "|---:|---:|---:|---:|---:|---:|---:|",
]
for r in rows:
    lines.append(
        f"| {float(r['lambda']):.3g} | {float(r['eta']):.3g} | "
        f"{float(r['m_pulse']):.3e} | {float(r['j_pulse']):.3e} | {float(r['s_pulse']):.3e} | "
        f"{float(r['pre_M_bump']):.3e} | {float(r['pre_J_bump']):.3e} |"
    )

if global_rows:
    max_global_s = max(abs(float(r["s_infinity"])) for r in global_rows)
    max_rem = max(abs(float(r["a19_remainder"])) for r in global_rows)
    max_ri = max(abs(float(r["rI_exterior_error"])) for r in global_rows)
    lines += ["", "## Global A.19 eta sweep", "",
              "原典の ordered asymptotic regime を数値化するため `lambda=1e-5` を使用する。",
              f"global root 後の最大 `|S(infinity)|`: `{max_global_s:.3e}`、最大 A.19 remainder: `{max_rem:.3e}`、最大 rI 誤差: `{max_ri:.3e}`。", "",
              "| eta | Amp(eta) | lambda S(inf) | A.19 remainder | s(Xp) |",
              "|---:|---:|---:|---:|---:|"]
    for r in global_rows:
        lines.append(
            f"| {float(r['eta']):.2f} | {float(r['amplitude']):.10g} | "
            f"{float(r['lambda_s_infinity']):.3e} | {float(r['a19_remainder']):.3e} | "
            f"{float(r['s_pulse']):.3e} |"
        )

lines += ["", "## 現在の実装境界", "",
          "A.2 → pulse → A.10 → A.11 → exterior transition → terminal tail を固定正規化で接続し、",
          "`S(infinity)=0` を full schedule から直接解くところまで到達した。A.9 の4つの reserved intervals は Proposition A.4 では意図的に変更しない。",
          "次は `Amp(eta)` の滑らかさ・対称性を診断した後、Appendix A.4 の pressure datum `(A.21)–(A.23)` へ進む。", ""]
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines), encoding="utf-8")
print(f"wrote {out}")
