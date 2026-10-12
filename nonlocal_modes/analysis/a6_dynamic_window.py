"""
a6_dynamic_window.py —— 动态嬗变与对称性破缺的可观测窗口 (论文 Fig.5 候选)
=========================================================================
只读取 E4b 的时间序列与汇总表 (results/E4_dynamic_transmutation/<tag>/e4b_*_timeseries.csv, e4b_summary.csv):
  (a) 嬗变核半径 r_s(xi): 稳定区 (三角形 m=2, U=25; h=0.1 与 h=0.05) 与失稳区 (U=15) 的 LG 发射, 虚线 = 定态 r_s*
  (b) C_3 破缺权重 (非 mu 权重) vs xi: 稳定区 LG 发射 (代数增长, 网格强迫)、定态发射与 1% 噪声 (不增长)、
      失稳区 (指数增长, 斜率 = E3 线性增长率)、圆盘 m=3 对照 (舍入种子 -> xi≈290 破裂)
  (c) 可观测窗口: 每个发射的嬗变距离 xi_t50 到破缺距离 xi_b (稳定区为传播终点 600)
运行:  python a6_dynamic_window.py
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import style  # noqa: E402

ROOT = os.path.join(io_utils.PROJECT_ROOT, "results", "E4_dynamic_transmutation")


def ts(tag, run):
    f = os.path.join(ROOT, tag, f"{run}_timeseries.csv")
    if not os.path.exists(f):
        return None
    rows = io_utils.read_csv(f)
    d = {k: np.array([float(r[k]) if r.get(k, "") not in ("", None) else np.nan for r in rows])
         for k in ("z", "nonmu_weight", "r_s_fit", "deviation")}
    return d


def summary(tag):
    f = os.path.join(ROOT, tag, "e4b_summary.csv")
    return {r["run"]: r for r in io_utils.read_csv(f)} if os.path.exists(f) else {}


def fnum(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(io_utils.PROJECT_ROOT, "results", "analysis", a.tag))
    S = summary(a.tag)
    S5 = summary(f"{a.tag}_h005")
    style.apply()
    fig, axs = style.plt.subplots(1, 3, figsize=(13.5, 4.0), gridspec_kw=dict(width_ratios=[1.0, 1.15, 1.0]))

    # (a) r_s(xi)
    ax = axs[0]
    series = [(a.tag, "e4b_poly3_m2_U25", "stable regime U=25 (h=0.1)", 0),
              (f"{a.tag}_h005", "e4b_poly3_m2_U25", "stable regime U=25 (h=0.05)", 1),
              (a.tag, "e4b_poly3_m2_U15", "unstable regime U=15 (h=0.1)", 2)]
    for tag, run, lab, i in series:
        d = ts(tag, run)
        if d is None:
            continue
        st = style.series_style(i)
        sel = d["z"] <= 120
        ax.plot(d["z"][sel], d["r_s_fit"][sel], "-", color=st["color"], lw=1.1, label=lab)
        r = (S5 if "h005" in tag else S).get(run)
        if r is not None and np.isfinite(fnum(r, "r_s_stationary")):
            ax.axhline(fnum(r, "r_s_stationary"), color=st["color"], lw=0.9, ls="--")
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("transmutation core radius r_s (harmonic fit)")
    ax.set_title("(a) LG m=2 into a triangle: the core forms within xi~2-5", loc="left")
    ax.text(0.98, 0.98, "dashed: stationary r_s*", transform=ax.transAxes, ha="right", va="top", fontsize=7,
            color=style.INK_SECONDARY)
    ax.set_ylim(-0.01, 0.42)
    ax.legend(loc="upper left", fontsize=7)

    # (b) 非 mu 权重
    ax = axs[1]
    cases = [(a.tag, "e4b_poly3_m2_U25", "stable, LG launch", 0),
             (a.tag, "e4b_poly3_m2_U25_stat", "stable, stationary launch", 1),
             (a.tag, "e4b_poly3_m2_U25_noise0.01", "stable, LG + 1% noise", 2),
             (a.tag, "e4b_poly3_m2_U15_noise1e-05", "unstable (U=15), noise 1e-5", 3),
             (a.tag, "e4b_disk_m3_U20", "disk m=3 control (round-off seed)", 4)]
    for tag, run, lab, i in cases:
        d = ts(tag, run)
        if d is None:
            continue
        st = style.series_style(i)
        ax.semilogy(d["z"], np.maximum(d["nonmu_weight"], 1e-32), "-", color=st["color"], lw=1.1, label=lab)
    for run, i in (("e4b_poly3_m2_U15_noise1e-05", 3), ("e4b_disk_m3_U20", 4)):
        r, d = S.get(run), ts(a.tag, run)
        if r is None or d is None:
            continue
        g = fnum(r, "linear_growth_E3_or_E5c")
        zb = fnum(r, "xi_b")
        if np.isfinite(g) and np.isfinite(zb):
            zz = np.linspace(max(zb - 120, 0), zb, 20)
            ax.semilogy(zz, 1e-3 * np.exp(2 * g * (zz - zb)), ":", color=style.INK, lw=1.2)
    ax.axhline(1e-3, color=style.INK_SECONDARY, lw=0.8)
    ax.text(600, 6e-4, "breaking threshold 1e-3", fontsize=7, color=style.INK_SECONDARY, va="top", ha="right")
    ax.text(0.98, 0.04, "dotted: exp(2 gamma xi) with the linear growth rate\nof the stationary state (E3)",
            transform=ax.transAxes, ha="right", fontsize=7, color=style.INK_SECONDARY)
    ax.set_ylim(1e-32, 3)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("weight outside the input irrep (C_n breaking)")
    ax.set_title("(b) symmetry breaking = seeded linear instability", loc="left")
    ax.legend(loc="center right", fontsize=6.8)

    # (c) 可观测窗口
    ax = axs[2]
    runs = [("e4b_poly3_m2_U25", "triangle U=25 (stable)"), ("e4b_poly3_m2_U30", "triangle U=30 (stable)"),
            ("e4b_poly3_m2_U15", "triangle U=15, grid seed"), ("e4b_poly3_m2_U15_noise1e-05", "triangle U=15, noise 1e-5"),
            ("e4b_heat3_m2_U15", "C3-heated disk, Theta=3"),
            ("e4b_heat2_m2_U15", "C2-heated disk, Theta=0.5 (twin)")]
    rows_out = []
    seen = set()
    for k, (run, lab) in enumerate(runs):
        r = S.get(run)
        if r is None:
            continue
        zt = fnum(r, "xi_t50")
        zb = fnum(r, "xi_b")
        zmax = fnum(r, "zmax")
        broken = np.isfinite(zb)
        end = zb if broken else zmax
        st = style.series_style(0 if not broken else 1)
        ax.plot([max(zt, 0.25), end], [k, k], "-", lw=6, color=st["color"], solid_capstyle="butt",
                label=("C_n symmetry kept to xi=600" if not broken else "C_n symmetry broken at xi_b")
                if (broken not in seen) else None)
        seen.add(broken)
        ax.plot(max(zt, 0.25), k, "o", color=style.INK, ms=4)
        if not broken:
            ax.annotate("", (end * 1.6, k), (end, k), arrowprops=dict(arrowstyle="->", color=st["color"]))
        rows_out.append(dict(run=run, xi_t50=zt, xi_b=zb, zmax=zmax))
        ax.text(0.3, k + 0.28, lab, fontsize=7, color=style.INK)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 1500)
    ax.set_yticks([])
    ax.set_ylim(-1.5, len(runs) + 1.0)
    ax.set_xlabel("propagation distance xi")
    ax.set_title("(c) observable window xi_t50 < xi < xi_b", loc="left")
    ax.text(0.98, 0.02, "dot: transmutation distance xi_t50\nbar end: breaking distance xi_b",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color=style.INK_SECONDARY)
    ax.legend(loc="upper left", fontsize=7)
    io_utils.write_csv(os.path.join(out, "observable_window.csv"), rows_out)
    style.save(fig, os.path.join(out, "dynamic_window.png"))
    print(f"written to {out}")


if __name__ == "__main__":
    main()
