"""
a3_e4_growth.py —— E4 长距离传播中 C_n 对称性破缺的增长率 vs E3 线性稳定性增长率
================================================================================
问题: E4 (xi=400) 中角赝动量 mu 的权重为什么大面积漂移? 是数值误差累积, 还是线性不稳定性把
网格误差/舍入误差这类“种子”放大? 后者的判据: 非 mu 分量的权重按 exp(2 gamma xi) 指数增长,
且 gamma 与 E3 中对应定态的增长率同量级 (E3 的 max Re(lambda) 包含保持 C_n 对称的模式, 所以是上界)。

做法 (只读已上传的 CSV):
  * 由 E4 的 <run>_irreps.csv 拟合非 mu 权重 W(xi) 在 1e-14 < W < 1e-2 段的指数增长: W ∝ exp(2 gamma xi);
    W 首次超过 1e-3 的位置记为破缺距离 xi_b;
  * 由 E2a 解族 (results/E2_transmutation_atlas/<tag>/<sample>_m<m>_summary.csv) 插值得到与输入能流 U
    对应的定态 b(U), 再由 E3 谱 (results/E3_stability/<tag>/<sample>_m<m>_spectrum.csv) 插值得到该 b 处的增长率。

输出 (results/analysis/<tag>/):  e4_growth_vs_e3.csv, e4_growth_vs_e3.png
运行:  python a3_e4_growth.py [--tag hq]
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import style  # noqa: E402


def fit_growth(z, w):
    ok = (w > 1e-14) & (w < 1e-2)
    zb = float(z[np.argmax(w > 1e-3)]) if (w > 1e-3).any() else float("nan")
    if ok.sum() < 5:
        return float("nan"), zb
    return float(np.polyfit(z[ok], np.log(w[ok]), 1)[0] / 2), zb


def interp_safe(x, xs, ys):
    o = np.argsort(xs)
    xs, ys = np.asarray(xs)[o], np.asarray(ys)[o]
    if x < xs[0] or x > xs[-1]:
        return float("nan")
    return float(np.interp(x, xs, ys))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    a = p.parse_args()
    root = io_utils.PROJECT_ROOT
    out = io_utils.ensure_dir(os.path.join(root, "results", "analysis", a.tag))
    rows = []
    for f in sorted(glob.glob(os.path.join(root, "results", "E4_dynamic_transmutation", a.tag, "*_irreps.csv"))):
        kind, m, U = re.findall(r"([a-z0-9]+)_m(\d+)_U([\d.]+)_irreps", os.path.basename(f))[0]
        m, U = int(m), float(U)
        tab = io_utils.read_csv(f)
        z = np.array([float(r["z"]) for r in tab])
        keys = [k for k in tab[0] if k.startswith("mu")]
        mu = max(keys, key=lambda k: float(tab[0][k]))
        w_other = np.array([sum(float(r[k]) for k in keys if k != mu) for r in tab])
        gam, zb = fit_growth(z, w_other)
        fam = os.path.join(root, "results", "E2_transmutation_atlas", a.tag, f"{kind}_m{m}_summary.csv")
        spec = os.path.join(root, "results", "E3_stability", a.tag, f"{kind}_m{m}_spectrum.csv")
        bU = g3 = float("nan")
        if os.path.exists(fam) and os.path.exists(spec):
            fr = io_utils.read_csv(fam)
            bU = interp_safe(U, [float(r["U"]) for r in fr], [float(r["b"]) for r in fr])
            sr = io_utils.read_csv(spec)
            if np.isfinite(bU):
                g3 = interp_safe(bU, [float(r["b"]) for r in sr], [float(r["max_growth"]) for r in sr])
        rows.append(dict(run=f"{kind}_m{m}_U{U:g}", sample=kind, m=m, U=U, b_of_U=bU, e3_growth_at_b=g3,
                         e4_symbreak_growth=gam, xi_symbreak=zb, max_nonmu_weight=float(w_other.max())))
        print(f"{kind:5s} m={m} U={U:4.0f}: b(U)={bU:5.2f}  E3 growth {g3:8.4f}  E4 symmetry-breaking growth "
              f"{gam:8.4f}  xi_b {zb:6.1f}  max non-mu weight {w_other.max():.1e}")
    io_utils.write_csv(os.path.join(out, "e4_growth_vs_e3.csv"), rows)

    style.apply()
    fig, ax = style.plt.subplots(figsize=(4.8, 4.0))
    for kind in style.SAMPLE_ORDER:
        sel = [r for r in rows if r["sample"] == kind and np.isfinite(r["e4_symbreak_growth"])
               and np.isfinite(r["e3_growth_at_b"])]
        if not sel:
            continue
        st = style.sample_style(kind)
        ax.plot([r["e3_growth_at_b"] for r in sel], [r["e4_symbreak_growth"] for r in sel], ls="none",
                color=st["color"], marker=st["marker"], ms=7, label=style.SAMPLE_LABEL[kind])
    lim = max(max(r["e3_growth_at_b"] for r in rows if np.isfinite(r["e3_growth_at_b"])), 0.12) * 1.05
    ax.plot([0, lim], [0, lim], "-", color=style.INK_SECONDARY, lw=1.0)
    ax.text(lim * 0.95, lim * 0.9, "y = x", ha="right", fontsize=7.5, color=style.INK_SECONDARY)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("E3: growth rate of the stationary state at b(U)")
    ax.set_ylabel("E4: growth rate of the C_n-breaking part")
    ax.set_title("C_n symmetry breaking in long propagation\nvs linear growth of the stationary state", loc="left")
    ax.legend(loc="upper left")
    style.save(fig, os.path.join(out, "e4_growth_vs_e3.png"))
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
