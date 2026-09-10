#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fv(row, key):
    return float(row[key])


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: generate_appendix_a4_report.py INPUT.csv OUTPUT.md")
    src, dst = sys.argv[1:]
    rows = load(src)
    if not rows:
        raise SystemExit("empty Appendix A.4 diagnostics")

    by_eta = {fv(r, "eta"): r for r in rows}
    even_err = 0.0
    odd_err = 0.0
    for eta, r in by_eta.items():
        if -eta in by_eta:
            rm = by_eta[-eta]
            even_err = max(even_err, abs(fv(r, "pi0_over_Pstar2") - fv(rm, "pi0_over_Pstar2")))
            odd_err = max(odd_err, abs(fv(r, "dpi0_deta_over_Pstar2") + fv(rm, "dpi0_deta_over_Pstar2")))

    min_margin = min(fv(r, "a22_margin") for r in rows)
    monotonic = [fv(r, "eta_times_dpi0_deta_over_Pstar2") for r in rows if abs(fv(r, "eta")) > 1e-14]
    min_sign = min(monotonic)
    center = min(rows, key=lambda r: abs(fv(r, "eta")))

    lines = [
        "# Appendix A.4 pressure datum 数値レポート",
        "",
        "原典の (A.21)–(A.23) に従い、A.2/A.3 の radial schedule から pressure datum を積分した軽量回帰。",
        "`Pi0/P*^2` を基本量にして、巨大な `P*` による不要な桁落ちを避けている。",
        "",
        "## パラメータ",
        "",
        f"- `lambda = {fv(center, 'lambda'):.8g}`",
        f"- `log(P*) = {fv(center, 'log_Pstar'):.8g}`",
        f"- `log(h) = {fv(center, 'log_h'):.8g}`",
        "- `Tf = 20`, `co = 0.05`",
        "",
        "この有限パラメータ計算は論文の存在証明そのものを数値証明するものではなく、式と schedule の実装回帰である。",
        "",
        "## 回帰結果",
        "",
        f"- `Pi0/P*^2` の偶対称性最大誤差: `{even_err:.3e}`",
        f"- `dPi0/deta/P*^2` の奇対称性最大誤差: `{odd_err:.3e}`",
        f"- (A.22) 内側 branch 境界の最小余裕: `{min_margin:.6e}`",
        f"- `eta * Pi0'(eta) / P*^2` の非零 eta での最小値: `{min_sign:.6e}`",
        f"- `eta=0` の `Pi0/P*^2`: `{fv(center, 'pi0_over_Pstar2'):.12e}`",
        "",
        "(A.22) の `-(5/2) P*^2 f(eta)^2` は full datum の等式ではなく、`y<=0` の reference inner branch だけから来る上界として検証している。",
        "A.11 の2つの angular E bump は total pressure increment を保存するため datum 積分からは省略し、A.11 区間そのものは保持した。",
        "",
        "## eta sweep",
        "",
        "| eta | Pi0/P*^2 | d_eta Pi0/P*^2 | A.22 margin | eta d_eta Pi0/P*^2 |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in rows:
        lines.append(
            f"| {fv(r,'eta'):.2f} | {fv(r,'pi0_over_Pstar2'):.9e} | "
            f"{fv(r,'dpi0_deta_over_Pstar2'):.9e} | {fv(r,'a22_margin'):.3e} | "
            f"{fv(r,'eta_times_dpi0_deta_over_Pstar2'):.3e} |"
        )
    lines += [
        "",
        "## (A.23) との接続",
        "",
        "inner branch `y<=0` では",
        "",
        "`Pi(y,eta)/P*^2 = Pi0(eta)/P*^2 + (5/2) f(eta)^2 exp(y/5)`",
        "",
        "を厳密形としてコード化した。`y -> -infinity` で `Pi0` に戻り、`d_y Pi = E^2/2 > 0` と整合する。",
        "",
        "次段階はこの `Pi0`, `Pi0_eta` を Appendix B の `Z*` に接続し、B.1–B.3 の analytic axis profile / moment matching を数値化する。",
        "",
    ]
    Path(dst).write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
