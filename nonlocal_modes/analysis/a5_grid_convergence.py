"""
a5_grid_convergence.py —— 三角形 m=2 失稳窗口上端的网格收敛 (论文 SI 图)
======================================================================
读取 E3b 四档网格 (h = 0.15 / 0.1 / 0.075 / 0.05; tag hq_coarse / hq / hq_h0075 / hq_h005) 的 e3b_summary.csv
与 a4_window_mechanism.py 写出的 hh_edge_check.csv (同一 b 点的增长率与由 Hamilton-Hopf 推出的耦合 g)。

三幅子图:
  (a) 窗口上端 b_c vs h^2 (误差棒 = 阈值括号), 线性外推到 h -> 0; 空心点 = 用该网格的耦合 g(h) 由 dw_disk(b) = 2g 预言的端点
  (b) 固定 b (2.2, 2.6, 2.8) 处的增长率 vs h^2 及线性外推
  (c) 推出的耦合 g vs h^2 (每档网格窗口内各点的中位数) —— 阈值随网格移动 = 粗网格低估了嬗变核的耦合
输出: results/analysis/<tag>/triangle_convergence.png、triangle_convergence.csv
运行:  python a5_grid_convergence.py    (先运行 a4_window_mechanism.py)
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import style  # noqa: E402

ROOT = os.path.join(io_utils.PROJECT_ROOT, "results")
GRIDS = [("hq_coarse", 0.15), ("hq", 0.1), ("hq_h0075", 0.075), ("hq_h005", 0.05)]


def upper_threshold(tag):
    f = os.path.join(ROOT, "E3_stability", tag, "e3b_summary.csv")
    if not os.path.exists(f):
        return None
    for r in io_utils.read_csv(f):
        if r["sample"] == "poly3" and r["m"] == "2" and float(r.get("theta0_deg", 0) or 0) == 0:
            thr = [t for t in r["thresholds"].split(";") if "unstable->stable" in t]
            if thr:
                lo, hi = thr[-1].split("(")[0].split("-")
                return float(lo), float(hi)
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    hh_file = os.path.join(out, "hh_edge_check.csv")
    if not os.path.exists(hh_file):
        sys.exit("run a4_window_mechanism.py first (hh_edge_check.csv)")
    hh = io_utils.read_csv(hh_file)
    rows = []
    for tag, h in GRIDS:
        thr = upper_threshold(tag)
        g_h = [float(r["g_inferred"]) for r in hh if abs(float(r["grid_h"]) - h) < 1e-9
               and r["source"] == "e3b_poly3_m2"]
        rows.append(dict(tag=tag, h=h, h2=h * h, bc_lo=thr[0] if thr else np.nan, bc_hi=thr[1] if thr else np.nan,
                         g_median=float(np.median(g_h)) if g_h else np.nan))
    h2 = np.array([r["h2"] for r in rows])
    bc = np.array([0.5 * (r["bc_lo"] + r["bc_hi"]) for r in rows])
    err = np.array([0.5 * (r["bc_hi"] - r["bc_lo"]) for r in rows])
    ok = np.isfinite(bc)
    fit_all = np.polyfit(h2[ok], bc[ok], 1)
    fit_fine = np.polyfit(h2[ok][-2:], bc[ok][-2:], 1)
    for r in rows:
        r.update(bc_extrap_all=fit_all[1], bc_extrap_two_finest=fit_fine[1])
    io_utils.write_csv(os.path.join(out, "triangle_convergence.csv"), rows)
    print(f"triangle m=2 upper edge: {', '.join(f'h={r['h']}: {r['bc_lo']:.4f}-{r['bc_hi']:.4f}' for r in rows)}")
    print(f"  extrapolation h->0: all four grids {fit_all[1]:.3f}, two finest {fit_fine[1]:.3f}")

    # Hamilton-Hopf 预言的端点: 每档网格用 g(h), 在圆盘频差 dw(b) (通道谱, b 到 3.4) 上找 dw = 2 g
    import a4_window_mechanism as a4
    ch = a4.load_channels()
    b1, w1 = a4.track(ch, -3, +1, 3.4, 0.886, max_jump=0.05)
    b0, w0 = a4.track(ch, 0, -1, 3.4, 0.848, max_jump=0.05)
    bvals = np.array([b for b in b1 if b >= 2.0 and np.min(np.abs(b0 - b)) < 1e-6])
    dwvals = np.array([np.interp(b, b1, w1) - np.interp(b, b0, w0) for b in bvals])

    style.apply()
    fig, axs = style.plt.subplots(1, 3, figsize=(12.5, 3.8))
    ax = axs[0]
    ax.errorbar(h2[ok], bc[ok], yerr=err[ok], fmt="o", color=style.PALETTE[0], mec=style.PALETTE[0], ms=6,
                capsize=3, label="E3b threshold (bracket)")
    xx = np.linspace(0, h2.max() * 1.05, 50)
    ax.plot(xx, np.polyval(fit_all, xx), "-", color=style.PALETTE[0], lw=1.0,
            label=f"linear in h^2: b_c(0) = {fit_all[1]:.2f}")
    ax.plot(xx, np.polyval(fit_fine, xx), "--", color=style.INK_SECONDARY, lw=1.0,
            label=f"two finest grids: b_c(0) = {fit_fine[1]:.2f}")
    for r in rows:
        if np.isfinite(r["g_median"]):
            # g 在窗口内略随 b 减小, 端点附近用该网格窗口内最靠右的点的 g
            gs = [float(x["g_inferred"]) for x in hh if abs(float(x["grid_h"]) - r["h"]) < 1e-9
                  and x["source"] == "e3b_poly3_m2"]
            bs = [float(x["b"]) for x in hh if abs(float(x["grid_h"]) - r["h"]) < 1e-9
                  and x["source"] == "e3b_poly3_m2"]
            g_edge = gs[int(np.argmax(bs))]
            if dwvals.min() < 2 * g_edge < dwvals.max():
                ax.plot(r["h2"], float(np.interp(2 * g_edge, dwvals, bvals)), "D", ms=7, color=style.PALETTE[1],
                        mfc="none", mec=style.PALETTE[1], mew=1.2,
                        label="Hamilton-Hopf edge, dw_disk = 2 g(h)" if r is rows[0] else None)
    ax.set_xlabel("grid spacing squared  h^2")
    ax.set_ylabel("upper edge of the window  b_c")
    ax.set_title("(a) triangle m=2: window edge vs grid", loc="left")
    ax.legend(loc="lower left", fontsize=7)

    ax = axs[1]
    for i, b0 in enumerate((2.2, 2.6, 2.8)):
        pts = sorted({float(x["grid_h"]): float(x["growth"]) for x in hh if abs(float(x["b"]) - b0) < 1e-6}.items())
        if len(pts) < 2:
            continue
        hs = np.array([p[0] ** 2 for p in pts])
        gs = np.array([p[1] for p in pts])
        st = style.series_style(i)
        ax.plot(hs, gs, st["marker"], ms=6, color=st["color"], mec=st["color"], label=f"b = {b0}")
        f2 = np.polyfit(hs[-2:], gs[-2:], 1)
        ax.plot([0, hs.max()], np.polyval(f2, [0, hs.max()]), "-", color=st["color"], lw=0.9)
        ax.plot(0, f2[1], st["marker"], ms=6, color=st["color"], mfc="none", mec=st["color"])
    ax.set_xlabel("grid spacing squared  h^2")
    ax.set_ylabel("growth rate of the kappa=0 mode")
    ax.set_ylim(bottom=0)
    ax.set_title("(b) growth at fixed b (open: extrapolated)", loc="left")
    ax.legend(loc="upper right", fontsize=7)

    ax = axs[2]
    gm = np.array([r["g_median"] for r in rows])
    okg = np.isfinite(gm)
    ax.plot(h2[okg], gm[okg], "o-", color=style.PALETTE[2], mec=style.PALETTE[2], ms=6)
    ax.set_xlabel("grid spacing squared  h^2")
    ax.set_ylabel("coupling g = sqrt(growth^2 + (dw/2)^2)")
    ax.set_title("(c) inferred C3 coupling of the L=-3 / L=0 pair", loc="left")
    style.save(fig, os.path.join(out, "triangle_convergence.png"))
    print(f"written to {out}")


if __name__ == "__main__":
    main()
