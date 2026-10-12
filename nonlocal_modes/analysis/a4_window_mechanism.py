"""
a4_window_mechanism.py —— 失稳窗口的机理图 (论文 Fig.4 候选): 圆盘 Bogoliubov 通道的 Krein 碰撞
================================================================================================
只读取已上传的 E3b / E5c 结果表:
  results/E3_stability/hq/e3b_disk_m2_channels.csv          圆盘 m=2 的通道分解谱 (b <= 1.35)
  results/E3_stability/hq_hib/e3b_disk_m2_channels.csv      同上, b = 1.4-3.4 (第三轮)
  results/E3_stability/<tag>/e3b_<样品>_m2_spectrum.csv      各样品 m=2 的失稳窗口 (增长率、失稳模频率)
  results/E5_thermal_reconfiguration/hq_smallT/e5c_n3_m1_T*_spectrum.csv   小 Theta 加热 (线性律)

四幅子图:
  (a) 低 b 段: 圆盘中 L=-2 (Krein +) 一支与 L=+2、+4、+6 (Krein -) 相交; 六边形 (伙伴 +4)、正方形 (伙伴 +2、+6) 的
      失稳模频率正好落在交叉处; 五边形的伙伴 L=+3 远在下方 -> 五边形窗口与圆盘相同。
  (b) 高 b 段: 三角形宽窗口的 kappa=0 对 L=-3 (K+) / L=0 (K-) 一直近简并; 三角形失稳模频率 = 两支平均值;
      Hamilton-Hopf: 由窗口内增长率推出耦合 g = sqrt(gamma^2 + (dw/2)^2), 窗口在 dw = 2g 处关闭。
  (c) 各样品 m=2 的增长率 vs b (窗口被拉宽的程度)。
  (d) 小 Theta 加热 (C_3, m=1): gamma/Theta vs b 塌缩到同一条半椭圆 (Hamilton-Hopf 线性律)。
另输出 hh_edge_check.csv: 三角形 m=2 各 b 点的 dw、gamma (各网格)、推出的 g 与预言的窗口端点。

运行:  python a4_window_mechanism.py           (输出到 results/analysis/hq/)
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

ROOT = os.path.join(io_utils.PROJECT_ROOT, "results")
# 通道的固定配色 (颜色跟随 L, 不随子图改变); Krein 符号用标记形状区分 (实心圆 +, 空心三角 -)
RHO_RECT = 12.638           # 20x40 矩形中心的保角半径 (数值解 Laplace 方程, 5 位收敛; 与正方形的 SC 公式交叉验证)
L_COLOR = {-2: style.PALETTE[0], 4: style.PALETTE[1], 2: style.PALETTE[2], 6: style.PALETTE[3],
           3: style.PALETTE[4], -3: style.PALETTE[5], 0: style.PALETTE[6], 1: style.PALETTE[7]}


def load_channels():
    rows = []
    for tag in ("hq", "hq_hib"):
        f = os.path.join(ROOT, "E3_stability", tag, "e3b_disk_m2_channels.csv")
        if os.path.exists(f):
            rows += [dict(b=float(r["b"]), im=float(r["lam_im"]), re=float(r["lam_re"]), L=int(r["L"]),
                          K=int(r["krein"])) for r in io_utils.read_csv(f)]
    return rows


def load_spectrum(tag, name):
    f = os.path.join(ROOT, "E3_stability", tag, f"{name}_spectrum.csv")
    if not os.path.exists(f):
        return []
    rows = [dict(b=float(r["b"]), g=float(r["growth"]), im=float(r["lam_im"]), unstable=r["unstable"] == "True",
                 L=r.get("L", ""), kappa=r.get("kappa", "")) for r in io_utils.read_csv(f)]
    return sorted(rows, key=lambda r: r["b"])


def branch(ch, L, K, im_lo, im_hi):
    """通道谱中 (L, Krein) 一支在频率带 [im_lo, im_hi] 内的点 (每个 b 取一个: 带内最接近带中心者)。"""
    out = {}
    for r in ch:
        if r["L"] == L and r["K"] == K and im_lo <= r["im"] <= im_hi:
            c = 0.5 * (im_lo + im_hi)
            if r["b"] not in out or abs(r["im"] - c) < abs(out[r["b"]] - c):
                out[r["b"]] = r["im"]
    bs = np.array(sorted(out))
    return bs, np.array([out[b] for b in bs])


def track(ch, L, K, b0, w0, max_jump=0.06):
    """沿 b 追踪 (L, Krein) 的一支: 从 (b0, w0) 出发向两侧, 每步取频率最接近上一点者 (跳变超过 max_jump 即停止)。"""
    pts = {}
    for r in ch:
        if r["L"] == L and r["K"] == K:
            pts.setdefault(round(r["b"], 4), []).append(r["im"])
    bs = sorted(pts)
    if not bs:
        return np.array([]), np.array([])
    i0 = int(np.argmin([abs(b - b0) for b in bs]))
    out = {}
    for direction in (+1, -1):
        w = w0
        i = i0
        while 0 <= i < len(bs):
            cand = pts[bs[i]]
            j = int(np.argmin([abs(c - w) for c in cand]))
            if abs(cand[j] - w) > max_jump:
                break
            out[bs[i]] = cand[j]
            w = cand[j]
            i += direction
    b = np.array(sorted(out))
    return b, np.array([out[x] for x in b])


def hh_check(ch, out):
    """三角形 m=2 的 kappa=0 窗口: 由圆盘两支的频差 dw 与各网格的增长率推出耦合 g, 并预言窗口端点 dw(b_e) = 2 g。"""
    b1, w1 = branch(ch, -3, +1, 0.45, 1.0)
    b0, w0 = branch(ch, 0, -1, 0.45, 1.0)
    common = sorted(set(np.round(b1, 4)) & set(np.round(b0, 4)))
    dw = {b: float(np.interp(b, b1, w1) - np.interp(b, b0, w0)) for b in common}
    wbar = {b: 0.5 * float(np.interp(b, b1, w1) + np.interp(b, b0, w0)) for b in common}
    grids = [("hq_coarse", 0.15), ("hq", 0.1), ("hq_h0075", 0.075), ("hq_h005", 0.05)]
    rows = []
    for tag, h in grids:
        for name in ("e3b_poly3_m2", "e3b_poly3_m2_rot15"):
            for r in load_spectrum(tag, name):
                if r["unstable"] and r["b"] >= 2.15 and r.get("kappa", "") in ("0", 0):
                    d = float(np.interp(r["b"], b1, w1) - np.interp(r["b"], b0, w0))
                    rows.append(dict(grid_h=h, source=name, b=r["b"], growth=r["g"], im_mode=r["im"],
                                     im_disk_mean=0.5 * float(np.interp(r["b"], b1, w1) + np.interp(r["b"], b0, w0)),
                                     dw_disk=d, g_inferred=float(np.sqrt(r["g"] ** 2 + (d / 2) ** 2))))
    if rows:
        g_fine = np.median([r["g_inferred"] for r in rows if r["grid_h"] <= 0.075])
        bb = np.array(common)
        dd = np.array([dw[b] for b in common])
        b_edge = float(np.interp(2 * g_fine, dd, bb)) if dd.min() < 2 * g_fine < dd.max() else float("nan")
        for r in rows:
            r["g_fine_median"] = g_fine
            r["b_edge_predicted"] = b_edge
        io_utils.write_csv(os.path.join(out, "hh_edge_check.csv"), rows)
        print(f"triangle m=2 kappa=0: coupling g (h<=0.075, median) = {g_fine:.4f}; predicted edge dw=2g at "
              f"b = {b_edge:.3f}")
        return dict(b=bb, dw=dd, wbar=np.array([wbar[b] for b in common]), g=g_fine, b_edge=b_edge, rows=rows)
    return None


def small_theta(out):
    data = []
    for f in sorted(glob.glob(os.path.join(ROOT, "E5_thermal_reconfiguration", "hq_smallT", "e5c_n3_m1_T*_spectrum.csv"))):
        T = float(re.findall(r"_T([\d.]+)_spectrum", f)[0])
        if T == 0:
            continue
        rows = sorted(io_utils.read_csv(f), key=lambda r: float(r["b"]))
        b = np.array([float(r["b"]) for r in rows])
        g = np.array([max(float(r["growth"]), 0.0) if r["unstable"] == "True" else 0.0 for r in rows])
        data.append((T, b, g))
    return sorted(data)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq", help="输出子目录 results/analysis/<tag>/")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    ch = load_channels()
    if not ch:
        sys.exit("no disk channel spectra found")
    hh = hh_check(ch, out)
    style.apply()
    fig, axs = style.plt.subplots(2, 2, figsize=(10.5, 7.6))

    # (a) 低 b 段的交叉: L=-2 (K+) 与 L=+2/+4/+6 (K-); 五边形的伙伴 L=+3 (K-)
    ax = axs[0, 0]
    for L, K, b0, w0 in ((-2, +1, 0.65, 0.450), (-2, -1, 0.65, 0.379), (2, -1, 0.65, 0.447), (4, -1, 0.65, 0.447),
                         (6, -1, 0.65, 0.441), (3, -1, 0.65, 0.210)):
        bs, ws = track(ch, L, K, b0, w0)
        sel = (bs >= 0.3) & (bs <= 1.36)
        col = L_COLOR[L]
        ax.plot(bs[sel], ws[sel], "o" if K > 0 else "^", ms=5, ls="-", lw=0.9, color=col,
                mfc=col if K > 0 else "none", mec=col, label=f"L={L:+d}, Krein {'+' if K > 0 else '-'}")
    for name, lab, dy in (("e3b_poly6_m2", "hexagon unstable mode", -10.0), ("e3b_poly4_m2", "square unstable mode", 6.0)):
        rows = [r for r in load_spectrum("hq", name) if r["unstable"] and r["b"] >= 0.5]
        if rows:
            ax.plot([r["b"] for r in rows], [r["im"] for r in rows], "-", color=style.INK, lw=5, alpha=0.18,
                    solid_capstyle="round")
            ax.annotate(lab, (rows[-1]["b"], rows[-1]["im"]), xytext=(8, -14 + dy), textcoords="offset points",
                        fontsize=7.5, color=style.INK_SECONDARY)
    ax.axvline(0.497, color=style.INK_SECONDARY, lw=0.8, ls=":")
    ax.text(0.505, 0.31, "disk window edge", fontsize=7, color=style.INK_SECONDARY)
    ax.set_xlim(0.3, 1.36)
    ax.set_ylim(0.15, 0.72)
    ax.text(0.56, 0.285, "pentagon partner L=+3 (Krein -): no crossing", fontsize=7, color=style.INK_SECONDARY)
    ax.set_xlabel("propagation constant b  (rho = 8)")
    ax.set_ylabel("frequency Im(lambda)")
    ax.set_title("(a) disk channels L and L+kn cross where polygon windows extend", loc="left")
    ax.legend(loc="upper left", fontsize=6.8, ncol=2)

    # (b) 三角形: L=-3 (K+) / L=0 (K-) 近简并 + 三角形失稳模频率
    ax = axs[0, 1]
    if hh is not None:
        b1, w1 = track(ch, -3, +1, 3.4, 0.886, max_jump=0.05)
        b0, w0 = track(ch, 0, -1, 3.4, 0.848, max_jump=0.05)
        ax.plot(b1, w1, "o", ms=4, ls="-", lw=0.9, color=L_COLOR[-3], mec=L_COLOR[-3], label="disk L=-3, Krein +")
        ax.plot(b0, w0, "^", ms=4.5, ls="-", lw=0.9, color=L_COLOR[0], mfc="none", mec=L_COLOR[0],
                label="disk L=0, Krein -")
        fine = [r for r in hh["rows"] if r["grid_h"] == 0.05]
        if fine:
            ax.plot([r["b"] for r in fine], [r["im_mode"] for r in fine], "s", ms=7, color=style.INK, mfc="none",
                    mec=style.INK, mew=1.2, ls="none", label="triangle unstable mode (h=0.05)")
        if np.isfinite(hh["b_edge"]):
            ax.axvline(hh["b_edge"], color=style.INK_SECONDARY, lw=0.9, ls="--")
            ax.text(hh["b_edge"] - 0.05, 0.36, f"predicted edge\n(dw = 2g, g={hh['g']:.4f})", fontsize=7,
                    color=style.INK_SECONDARY, ha="right")
    ax.set_xlabel("propagation constant b  (rho = 8)")
    ax.set_ylabel("frequency Im(lambda)")
    ax.set_title("(b) triangle: near-degenerate kappa=0 pair (L=-3 / L=0)", loc="left")
    ax.legend(loc="upper left", fontsize=7)

    # (c) 各样品 m=2 的增长率
    ax = axs[1, 0]
    for kind in ("disk", "poly5", "poly6", "poly4", "poly3", "rect"):
        rows = load_spectrum("hq", f"e3b_{kind}_m2")
        if not rows:
            continue
        st = style.sample_style(kind)
        f = (RHO_RECT / 8.0) ** 2 if kind == "rect" else 1.0     # 标度不变性: b、增长率都乘 (rho/8)^2
        ax.plot([f * r["b"] for r in rows], [f * max(r["g"], 0) for r in rows], "-", color=st["color"],
                marker=st["marker"], ms=3.5, mec=st["color"],
                label=style.SAMPLE_LABEL[kind] + (" 20x40, rescaled to rho=8" if kind == "rect" else ""))
    ax.set_xlim(0, 3.1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("propagation constant b")
    ax.set_ylabel("growth rate max Re(lambda)")
    ax.set_title("(c) m=2: the intrinsic window and its extensions", loc="left")
    ax.legend(loc="upper right", fontsize=7)

    # (d) 小 Theta 线性律: gamma/Theta 塌缩
    ax = axs[1, 1]
    data = small_theta(out)
    for i, (T, b, g) in enumerate(data):
        st = style.series_style(i)
        sel = g > 0
        ax.plot(b[sel], g[sel] / T, "-", color=st["color"], marker=st["marker"], ms=3.5, lw=1.2, mec=st["color"],
                label=f"Theta = {T:g}  (max/Theta = {g.max() / T:.4f})")
    ax.set_xlabel("propagation constant b  (R = 10)")
    ax.set_ylabel("growth rate / Theta")
    ax.set_title("(d) C3 heating, m=1: growth proportional to the coupling", loc="left")
    ax.legend(loc="lower center", fontsize=7)
    style.save(fig, os.path.join(out, "window_mechanism.png"))
    print(f"written to {out}")


if __name__ == "__main__":
    main()
