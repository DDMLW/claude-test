"""
a8_heated_window.py —— C_3 加热 m=2 高 b 失稳窗口 = 圆盘 κ=0 通道对 L=-3 / L=0 的 Hamilton-Hopf 碰撞, 耦合 ∝ Theta (论文 SI)
=====================================================================================================================
读取:
  results/E5_thermal_reconfiguration/hq_hib/e5c_n3_m2_T*_spectrum.csv   C_3 加热 m=2, b=2.9-6.0 (R=10, h=0.1)
  results/E5_thermal_reconfiguration/hq/e5c_n3_m2_T*_spectrum.csv       同上, b<=3.0 (第二轮)
  圆盘 (R=10, 与 E5 同一网格) m=2 的通道谱: --channels 给出的 csv (默认 results/E3_stability/hq_disk10/e3b_disk_m2_channels.csv,
      第四轮正式数据; 之前可用云端 scratch 计算的同格式文件)
二模 Hamilton-Hopf: gamma^2 = g^2 - (delta/2)^2, delta(b) = 圆盘 L=-3 (K+) 与 L=0 (K-) 两支的频差, 失稳模频率 = 两支平均 (+O(Theta^2));
若耦合 g = c(b) Theta, 则 c = sqrt(gamma^2 + delta^2/4) / Theta 对不同 Theta 应塌缩到同一条 c(b), 窗口上端满足 delta(b_up) = 2 c(b_up) Theta。
三幅子图:
  (a) 圆盘两支频率 vs b 与各 Theta 的失稳模频率
  (b) 推出的 c = sqrt(gamma^2 + delta^2/4)/Theta vs b (各 Theta 是否塌缩)
  (c) 预言的上端 b_up(Theta) (delta(b) = 2 c Theta) 对 E5c 实测
输出: results/analysis/<tag>/heated_window.png、heated_window.csv
运行:  python a8_heated_window.py [--channels <通道谱 csv>]
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


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def disk_pair(fname):
    """圆盘通道谱中 Im 0.55-1.3 段的 L=-3 (K+) 与 L=0 (K-) 两支 (每个 b 取该段最低的一个)。返回 {b: (w_m3, w_0)}。"""
    out = {}
    rows = io_utils.read_csv(fname)
    for b in sorted({round(fnum(r["b"]), 6) for r in rows}):
        sel = [r for r in rows if abs(fnum(r["b"]) - b) < 1e-9 and 0.55 <= fnum(r["lam_im"]) <= 1.3]
        wm3 = sorted(fnum(r["lam_im"]) for r in sel if r["L"] == "-3" and r["krein"] == "1")
        w0 = sorted(fnum(r["lam_im"]) for r in sel if r["L"] == "0" and r["krein"] == "-1")
        if wm3 and w0:
            out[b] = (wm3[0], w0[0])
    return out


def heated_rows():
    """各 Theta 的 (b, growth, Im) —— hq_hib 优先, 补充 hq 中 b<2.9 的点。"""
    data = {}
    for tag in ("hq", "hq_hib"):
        for f in glob.glob(os.path.join(ROOT, "E5_thermal_reconfiguration", tag, "e5c_n3_m2_T*_spectrum.csv")):
            T = float(re.search(r"_T([0-9.]+)_spectrum", f).group(1))
            if T == 0:
                continue
            for r in io_utils.read_csv(f):
                b = round(fnum(r["b"]), 6)
                if tag == "hq" and b >= 2.9 - 1e-9:
                    continue
                data.setdefault(T, {})[b] = (fnum(r["growth"]), fnum(r["lam_im"]), r.get("kappa", ""))
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    p.add_argument("--channels", default=os.path.join(ROOT, "E3_stability", "hq_disk10", "e3b_disk_m2_channels.csv"))
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    if not os.path.exists(a.channels):
        sys.exit(f"channel spectrum not found: {a.channels}")
    pair = disk_pair(a.channels)
    bd = np.array(sorted(pair))
    wm3 = np.array([pair[b][0] for b in bd])
    w0 = np.array([pair[b][1] for b in bd])
    delta = wm3 - w0
    data = heated_rows()
    Ts = sorted(data)
    rows = []
    for T in Ts:
        for b, (g, im, kap) in sorted(data[T].items()):
            if not (bd.min() <= b <= bd.max()) or g < 1e-5 or kap not in ("0", "0.0"):
                continue
            d = float(np.interp(b, bd, delta))
            wbar = float(np.interp(b, bd, 0.5 * (wm3 + w0)))
            rows.append(dict(Theta=T, b=b, growth=g, im_mode=im, disk_mean=wbar, delta_disk=d,
                             c_inferred=np.sqrt(g * g + d * d / 4) / T))
    # 每个 Theta 的实测上端 (E5c) 与预言的上端: delta(b) = 2 c Theta, c 取各 Theta 推出值的中位数 (随 b 插值)
    edges = {}
    for r in io_utils.read_csv(os.path.join(ROOT, "E5_thermal_reconfiguration", "hq_hib", "e5c_upper_edges.csv")):
        edges[fnum(r["Theta"])] = (fnum(r["b_up"]), r["b_up_is_bound"] == "True")
    bc = np.array(sorted({r["b"] for r in rows}))
    c_med = np.array([np.median([r["c_inferred"] for r in rows if abs(r["b"] - b) < 1e-9]) for b in bc])
    pred = {}
    for T in Ts:
        bb = np.linspace(bc.min(), bd.max(), 400)
        f = np.interp(bb, bd, delta) - 2 * T * np.interp(bb, bc, c_med)
        k = np.where((f[:-1] < 0) & (f[1:] >= 0))[0]
        pred[T] = float(bb[k[-1]] - f[k[-1]] * (bb[k[-1] + 1] - bb[k[-1]]) / (f[k[-1] + 1] - f[k[-1]])) if len(k) else np.nan
    for r in rows:
        r["b_up_measured"] = edges.get(r["Theta"], (np.nan, False))[0]
        r["b_up_predicted"] = pred.get(r["Theta"], np.nan)
    io_utils.write_csv(os.path.join(out, "heated_window.csv"), rows)
    print("C_3 heated m=2, high-b window (kappa=0 pair L=-3/L=0):")
    for T in Ts:
        sel = [r for r in rows if r["Theta"] == T]
        if sel:
            print(f"  Theta={T:g}: c = {np.median([r['c_inferred'] for r in sel]):.4f} "
                  f"(range {min(r['c_inferred'] for r in sel):.4f}-{max(r['c_inferred'] for r in sel):.4f}); "
                  f"max |Im - disk mean| = {max(abs(r['im_mode'] - r['disk_mean']) for r in sel):.4f}; "
                  f"b_up measured {edges.get(T, (np.nan,))[0]:.3f}, predicted {pred[T]:.3f}")

    style.apply()
    fig, axs = style.plt.subplots(1, 3, figsize=(13.0, 3.9))
    Ts = [T for T in Ts if any(r["Theta"] == T for r in rows)]       # 颜色只分给有数据的 Theta (固定顺序)
    ax = axs[0]
    # 圆盘两支是参照线: 用中性墨色 (实心圆 = Krein +, 空心三角 = Krein -), 不占用数据系列的颜色
    ax.plot(bd, wm3, "o-", color=style.INK, mec=style.INK, ms=4, lw=1.0, label="disk L=-3 (Krein +)")
    ax.plot(bd, w0, "^-", color=style.INK_SECONDARY, mfc="none", mec=style.INK_SECONDARY, ms=5, lw=1.0,
            label="disk L=0 (Krein -)")
    for i, T in enumerate(Ts):
        sel = [r for r in rows if r["Theta"] == T]
        if sel:
            st = style.series_style(i)
            ax.plot([r["b"] for r in sel], [r["im_mode"] for r in sel], st["marker"], ms=3.5, color=st["color"],
                    mec=st["color"], ls="none", label=f"unstable mode, Theta={T:g}")
    ax.set_xlabel("propagation constant b (R=10)")
    ax.set_ylabel("frequency Im(lambda)")
    ax.set_title("(a) heated-state mode sits between the disk pair", loc="left")
    ax.legend(loc="upper left", fontsize=6.5)

    ax = axs[1]
    for i, T in enumerate(Ts):
        sel = [r for r in rows if r["Theta"] == T]
        if sel:
            st = style.series_style(i)
            ax.plot([r["b"] for r in sel], [r["c_inferred"] for r in sel], st["marker"] + "-", ms=4, lw=0.9,
                    color=st["color"], mec=st["color"], label=f"Theta={T:g}")
    ax.set_xlabel("propagation constant b (R=10)")
    ax.set_ylabel("c = sqrt(gamma^2 + delta^2/4) / Theta")
    ax.set_ylim(bottom=0)
    ax.set_title("(b) coupling per unit heating collapses", loc="left")
    ax.legend(loc="lower left", fontsize=7)

    ax = axs[2]
    Tm = [T for T in Ts if T in edges and np.isfinite(pred.get(T, np.nan))]
    ax.plot(Tm, [edges[T][0] for T in Tm], "o", ms=7, color=style.PALETTE[0], mec=style.PALETTE[0],
            label="E5c upper edge b_up")
    ax.plot(Tm, [pred[T] for T in Tm], "D", ms=7, color=style.PALETTE[1], mfc="none", mec=style.PALETTE[1],
            label="Hamilton-Hopf: delta_disk(b) = 2 c Theta")
    ax.set_xlabel("heating amplitude Theta")
    ax.set_ylabel("upper edge of the window b_up")
    ax.set_title("(c) window edge set by the disk detuning", loc="left")
    ax.legend(loc="upper left", fontsize=7)
    style.save(fig, os.path.join(out, "heated_window.png"))
    print(f"written to {out}")


if __name__ == "__main__":
    main()
