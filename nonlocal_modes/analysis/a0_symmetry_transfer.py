"""
a0_symmetry_transfer.py —— 对称性传递定律与保角半径普适性 (论文 Fig.1)
=====================================================================
只读取 E1 的结果表 (results/E1_symmetry_transfer/<tag>/):
  e1a_index_anisotropy.csv   C_n 各向异性分量 dn_n(r) (part=power_law / universality / linearity / rectangle)
  e1b_all.csv                rho=8 的多边形与圆盘中 m=0、1 解族的 U(b)
四幅子图:
  (a) 塌缩: dn_n * pi n (n+1) / U  vs  r/rho (双对数), 各 n、各 rho 的点落在 (r/rho)^n 上 —— 定律 δn = U/(πn(n+1)) (r/ρ)^n cos nθ
  (b) 数值/理论 - 1: 不同光束形状 (高斯 w=1-3、环形) 与不同能流 U=2-40 —— 系数与光束形状无关、对 U 线性
  (c) rho=8 的三角形...六边形与圆盘 m=0 的 U(b) 重合 (保角半径普适性)
  (d) 相对圆盘的 |ΔU|/U vs b (m=0 实线、m=1 虚线)
输出: results/analysis/<tag>/symmetry_transfer.png、symmetry_transfer_summary.csv
运行:  python a0_symmetry_transfer.py
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import style  # noqa: E402

ROOT = os.path.join(io_utils.PROJECT_ROOT, "results")
POLY = ["poly3", "poly4", "poly5", "poly6"]


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    a = p.parse_args()
    src = os.path.join(ROOT, "E1_symmetry_transfer", a.tag)
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    e1a = io_utils.read_csv(os.path.join(src, "e1a_index_anisotropy.csv"))
    e1b = io_utils.read_csv(os.path.join(src, "e1b_all.csv"))
    summary = {}
    style.apply()
    fig, axs = style.plt.subplots(1, 4, figsize=(15.0, 3.7))

    # (a) 塌缩
    ax = axs[0]
    pl = [r for r in e1a if r["part"] == "power_law"]
    for i, n in enumerate((3, 4, 5, 6)):
        sel = [r for r in pl if int(r["n"]) == n]
        x = np.array([fnum(r["r"]) / fnum(r["rho"]) for r in sel])
        y = np.array([fnum(r["dn_n"]) * np.pi * n * (n + 1) / fnum(r["U"]) for r in sel])
        st = style.series_style(i)
        ok = y > 0
        ax.loglog(x[ok], y[ok], st["marker"], ms=4, color=st["color"], mec=st["color"], ls="none", label=f"n={n}")
        xx = np.array([x[ok].min(), x[ok].max()])
        ax.loglog(xx, xx ** n, "-", color=st["color"], lw=0.9)
        rat = np.array([fnum(r["ratio"]) for r in sel])
        summary[f"n{n}_ratio_min"], summary[f"n{n}_ratio_max"] = float(rat.min()), float(rat.max())
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.06, 0.1, 0.2, 0.4, 0.6]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("r / rho")
    ax.set_ylabel("dn_n * pi n (n+1) / U")
    ax.set_title("(a) transfer law collapses onto (r/rho)^n", loc="left")
    ax.text(0.98, 0.03, "lines: (r/rho)^n, no fit", transform=ax.transAxes, ha="right", fontsize=7,
            color=style.INK_SECONDARY)
    ax.legend(loc="upper left", fontsize=7)

    # (b) 普适性与线性
    ax = axs[1]
    uni = [r for r in e1a if r["part"] == "universality"]
    lin = [r for r in e1a if r["part"] == "linearity"]
    ax.plot(range(len(uni)), [fnum(r["ratio"]) - 1 for r in uni], "o", ms=6, color=style.PALETTE[0],
            mec=style.PALETTE[0], label="beam shape (U=10)")
    ax.plot(np.arange(len(lin)) + len(uni) + 1, [fnum(r["ratio"]) - 1 for r in lin], "s", ms=6,
            color=style.PALETTE[1], mec=style.PALETTE[1], label="power U (Gaussian w=1.5)")
    labels = [r["beam"].replace("gauss ", "G ").replace("ring ", "ring ") for r in uni] + [""] + \
             [f"U={fnum(r['U']):g}" for r in lin]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=6.5)
    ax.axhline(0, color=style.INK_SECONDARY, lw=0.8)
    ax.set_ylim(-5e-5, 5e-5)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.set_ylabel("numerical / theory - 1   (n=3, rho=10)")
    ax.set_title("(b) coefficient independent of beam shape", loc="left")
    ax.legend(loc="upper left", fontsize=7)
    summary["universality_max_dev"] = float(max(abs(fnum(r["ratio"]) - 1) for r in uni + lin))

    # (c) U(b) 重合
    ax = axs[2]
    fam = {}
    for r in e1b:
        fam.setdefault((r["sample"], int(r["m"])), []).append((fnum(r["b"]), fnum(r["U"])))
    # 圆盘先画 (粗、在下), 多边形用细虚线叠在上面: 五条曲线在图上完全重合, 差别见 (d)
    for kind in ["disk"] + POLY:
        st = style.sample_style(kind)
        d = np.array(sorted(fam[(kind, 0)]))
        ax.plot(d[:, 0], d[:, 1], "-" if kind == "disk" else "--", color=st["color"],
                lw=3.0 if kind == "disk" else 1.0, label=style.SAMPLE_LABEL[kind])
    ax.text(0.98, 0.03, "all five curves coincide\n(max deviation 3.7e-4, panel d)", transform=ax.transAxes,
            ha="right", fontsize=7, color=style.INK_SECONDARY)
    ax.set_xlabel("propagation constant b")
    ax.set_ylabel("power U   (m=0, rho=8)")
    ax.set_title("(c) same conformal radius, same U(b)", loc="left")
    ax.legend(loc="upper left", fontsize=7)

    # (d) 相对圆盘的差
    ax = axs[3]
    for m, ls in ((0, "-"), (1, "--")):
        ref = dict(fam[("disk", m)])
        for kind in POLY:
            st = style.sample_style(kind)
            d = sorted(fam[(kind, m)])
            bs = np.array([b for b, _ in d if b in ref])
            dev = np.array([abs(u - ref[b]) / ref[b] for b, u in d if b in ref])
            ax.semilogy(bs, dev, ls, color=st["color"], lw=1.0,
                        label=f"{style.SAMPLE_LABEL[kind].split(' ')[0]}" if m == 0 else None)
            sel = bs >= (0.5 if m == 1 else 0.0)
            summary[f"{kind}_m{m}_max_rel_dev"] = float(dev[sel].max()) if sel.any() else float("nan")
    ax.set_xlabel("propagation constant b")
    ax.set_ylabel("|U - U_disk| / U_disk")
    ax.set_title("(d) |U - U_disk| / U_disk", loc="left")
    ax.text(0.98, 0.97, "solid: m=0   dashed: m=1", transform=ax.transAxes, ha="right", va="top", fontsize=7,
            color=style.INK_SECONDARY)
    ax.legend(loc="lower left", fontsize=7, ncol=2)
    io_utils.write_csv(os.path.join(out, "symmetry_transfer_summary.csv"), [summary])
    style.save(fig, os.path.join(out, "symmetry_transfer.png"))
    for k, v in summary.items():
        print(f"  {k}: {v:.3e}")
    print(f"written to {out}")


if __name__ == "__main__":
    main()
