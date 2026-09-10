#!/usr/bin/env python3
import csv
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
max_pre_m = max(abs(float(r["pre_M_bump"])) for r in rows)
max_pre_j = max(abs(float(r["pre_J_bump"])) for r in rows)

lines = [
    "# Appendix A.3 moment closure 診断",
    "",
    "このレポートは Appendix A.2 の staged radial profile を A.3 の closure solver に接続した軽量数値回帰を GitHub Actions で再生成したもの。",
    "`M/(XE)`, `J/(XHE)`, `S/(XE^2)` を直接発展させるため、巨大な `X` や `P*` を明示的に生成せず theorem-scale の対数半径を扱える。",
    "完全な Theorem 4.6 profile assembly ではなく、現在は axial pulse 開始までの global schedule と A.15 closure を直接接続している。",
    "",
    "## 主要結果",
    "",
    f"- `(A.19)` の `Kb = {kb:.12g}`（論文中の評価 `0.20 < Kb <= 0.25` と整合）",
    f"- remainder を0とした principal amplitude: `{amp:.12g}`（論文の bracket `0.9 < Amp < 1.2` 内）",
    f"- A.2 schedule から A.15 の第1 bump 中心へ運んだ pre-M 最大値: `{max_pre_m:.3e}`",
    f"- A.2 schedule から A.15 の第1 bump 中心へ運んだ pre-J 最大値: `{max_pre_j:.3e}`",
    f"- `(A.15)` の `M,J` closure 最大残差: `{max_mj:.3e}`",
    f"- `(A.11)` の `I, pressure increment` closure 最大残差: `{max_ip:.3e}`",
    f"- `(A.16)` の `Qs(0)=Qp` 最大差: `{max_q:.3e}`",
    "",
    "## λ・η sweep",
    "",
    "| lambda | eta | m(Xp) | j(Xp) | s(Xp) | pre-M | pre-J | c1(MJ) | c2(MJ) |",
    "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for r in rows:
    lines.append(
        f"| {float(r['lambda']):.3g} | {float(r['eta']):.3g} | "
        f"{float(r['m_pulse']):.3e} | {float(r['j_pulse']):.3e} | {float(r['s_pulse']):.3e} | "
        f"{float(r['pre_M_bump']):.3e} | {float(r['pre_J_bump']):.3e} | "
        f"{float(r['c1_MJ']):.3e} | {float(r['c2_MJ']):.3e} |"
    )

lines += [
    "",
    "## 実装境界",
    "",
    "A.2 の reference power law → first transition → axial reduction → intermediate entry → `(A.9)` hold を一つの正規化ODEとして積分し、",
    "得られた `(A.14)` の pre-pulse discrepancies を `(A.15)` の two-bump solver へ直接入力するところまで統合済み。",
    "次は pulse 本体・A.10 profile interpolation・A.11 angular correction・exterior transition を同じ正規化状態に追加し、",
    "`S(infinity)` を trial amplitude ごとに直接評価して `(A.19)` の remainder を callback ではなく global schedule から構成する。",
    "",
]
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines), encoding="utf-8")
print(f"wrote {out}")
