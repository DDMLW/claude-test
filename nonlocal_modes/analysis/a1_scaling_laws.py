"""
a1_scaling_laws.py —— 标度律分析: r_s 随保角半径 rho (E2b) 和加热幅度 Theta (E5a) 的幂律指数
=============================================================================================
只读取已上传的结果表 (results/E2_transmutation_atlas/<tag>/scaling_n*_m*.csv,
results/E5_thermal_reconfiguration/<tag>/heating_n*_m*.csv), 不需要原始数据, 在任何电脑上几秒钟跑完。

理论 (docs/theory.md 第 3 节): 弱各向异性极限 r_s ∝ eps^s, s 由牛顿多边形给出;
固定 U 时 eps ∝ rho^(-n) -> r_s ∝ rho^(-n s) (E2b); eps ∝ Theta -> r_s ∝ Theta^s (E5a)。
这是渐近律, 只在卫星环远小于亮环 (r_s / ring_radius -> 0) 时成立, 所以除全局拟合外还给出:
  * 相邻两点的局部斜率, 以及最靠近渐近区 (r_s/ring 最小) 的局部斜率 slope_asym;
  * 渐近窗口 (r_s/ring <= 0.2, 至少 3 点) 内的拟合斜率 slope_window。

输出 (results/analysis/<tag>/):
  scaling_summary.csv        每种情形: 理论斜率、全局拟合、窗口拟合、渐近局部斜率、点数、r_s/ring 范围
  scaling_local_slopes.csv   每对相邻点的局部斜率
  scaling_local_slopes.png   局部斜率 / 理论斜率 vs r_s/ring (渐近区应趋于 1)

运行:  python a1_scaling_laws.py [--tag hq]
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils, theory  # noqa: E402

import style  # noqa: E402

WINDOW = 0.2       # 渐近窗口: r_s / ring_radius <= WINDOW


def load_cases(tag):
    root = io_utils.PROJECT_ROOT
    cases = []
    for f in sorted(glob.glob(os.path.join(root, "results", "E2_transmutation_atlas", tag, "scaling_n*_m*.csv"))):
        n, m = map(int, re.findall(r"scaling_n(\d+)_m(\d+)", f)[0])
        rows = io_utils.read_csv(f)
        x = np.array([float(r["rho"]) for r in rows])
        cases.append(dict(source="E2b", n=n, m=m, xname="rho", x=x, rows=rows))
    for f in sorted(glob.glob(os.path.join(root, "results", "E5_thermal_reconfiguration", tag, "heating_n*_m*.csv"))):
        n, m = map(int, re.findall(r"heating_n(\d+)_m(\d+)", f)[0])
        rows = [r for r in io_utils.read_csv(f) if float(r["Theta"]) > 0 and r.get("hybrid", "False") != "True"]
        x = np.array([float(r["Theta"]) for r in rows])
        cases.append(dict(source="E5a", n=n, m=m, xname="Theta", x=x, rows=rows))
    for c in cases:
        rows = c["rows"]
        c["rs"] = np.array([float(r["r_s_fit"]) for r in rows])
        c["ring"] = np.array([float(r["ring_radius"]) for r in rows])
        s = theory.selection_rule(c["m"], c["n"])["rings"][0][2]
        c["s"] = s
        c["slope_theory"] = -c["n"] * s if c["source"] == "E2b" else s
        ok = np.isfinite(c["rs"]) & (c["rs"] > 0)
        o = np.argsort(c["x"][ok])
        c["x"], c["rs"], c["ring"] = c["x"][ok][o], c["rs"][ok][o], c["ring"][ok][o]
    return cases


def analyse(c):
    x, rs, ring = c["x"], c["rs"], c["ring"]
    ratio = rs / ring
    loc = np.diff(np.log(rs)) / np.diff(np.log(x))
    mid = np.sqrt(ratio[1:] * ratio[:-1])
    glob_fit = np.polyfit(np.log(x), np.log(rs), 1)[0]
    i_asym = int(np.argmin(mid))
    win = ratio <= WINDOW
    if win.sum() < 3:
        win = np.zeros_like(win)
        win[np.argsort(ratio)[:3]] = True
    win_fit = np.polyfit(np.log(x[win]), np.log(rs[win]), 1)[0]
    return dict(local=loc, mid=mid, slope_global=glob_fit, slope_asym=loc[i_asym],
                ratio_asym=mid[i_asym], slope_window=win_fit, n_window=int(win.sum()),
                ratio_min=ratio.min(), ratio_max=ratio.max())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(io_utils.PROJECT_ROOT, "results", "analysis", a.tag))
    cases = load_cases(a.tag)
    if not cases:
        sys.exit(f"no scaling/heating tables found for tag {a.tag}")
    summary, local_rows = [], []
    for c in cases:
        r = analyse(c)
        c.update(r)
        name = f"{c['source']} (n,m)=({c['n']},{c['m']})"
        summary.append(dict(source=c["source"], n=c["n"], m=c["m"], variable=c["xname"], s_theory=round(c["s"], 4),
                            slope_theory=round(c["slope_theory"], 4), slope_global=r["slope_global"],
                            slope_window=r["slope_window"], n_window=r["n_window"], slope_asym=r["slope_asym"],
                            ratio_at_asym=r["ratio_asym"], rs_over_ring_min=r["ratio_min"],
                            rs_over_ring_max=r["ratio_max"], n_points=len(c["x"])))
        for xa, xb, sl, md in zip(c["x"][:-1], c["x"][1:], r["local"], r["mid"]):
            local_rows.append(dict(source=c["source"], n=c["n"], m=c["m"], x_from=xa, x_to=xb, rs_over_ring=md,
                                   local_slope=sl, local_over_theory=sl / c["slope_theory"]))
        print(f"{name:22s} theory {c['slope_theory']:+.3f} | global {r['slope_global']:+.3f} | "
              f"window(r_s/ring<={WINDOW}, {r['n_window']} pts) {r['slope_window']:+.3f} | "
              f"asymptotic local {r['slope_asym']:+.3f} at r_s/ring={r['ratio_asym']:.3f}")
    io_utils.write_csv(os.path.join(out, "scaling_summary.csv"), summary)
    io_utils.write_csv(os.path.join(out, "scaling_local_slopes.csv"), local_rows)

    style.apply()
    fig, axs = style.plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
    for ax, src, title in ((axs[0], "E2b", "(a) shape anisotropy: r_s vs conformal radius rho"),
                           (axs[1], "E5a", "(b) boundary heating: r_s vs Theta")):
        sel = [c for c in cases if c["source"] == src]
        for i, c in enumerate(sel):
            st = style.series_style(i)
            ax.semilogx(c["mid"], c["local"] / c["slope_theory"], "-", color=st["color"], marker=st["marker"],
                        label=f"(n,m)=({c['n']},{c['m']}): s={c['s']:.3g}")
        ax.axhline(1.0, color=style.INK_SECONDARY, lw=1.0)
        xmin = min(c["mid"].min() for c in sel) * 0.7
        ax.axvspan(xmin, WINDOW, color=style.GRID, alpha=0.5, lw=0)
        ax.text(WINDOW * 0.9, 1.3, "asymptotic window", ha="right", va="top", fontsize=7.5,
                color=style.INK_SECONDARY)
        ax.set_xlim(left=xmin)
        ax.set_xlabel("satellite radius / ring radius  (r_s / r_ring)")
        ax.set_title(title, loc="left")
        ax.legend(loc="lower left")
    axs[0].set_ylabel("local exponent / theory exponent")
    axs[0].set_ylim(0.0, 1.35)
    style.save(fig, os.path.join(out, "scaling_local_slopes.png"))
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
