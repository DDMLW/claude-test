"""
a2_stability_map.py —— 稳定性图谱: 各样品、各拓扑荷的不稳定区间与阈值
====================================================================
只读取已上传的 E3 结果 (results/E3_stability/<tag>/<sample>_m<m>_spectrum.csv), 几秒钟跑完。

对每个 (样品, m):
  * 不稳定 b 区间 (max Re lambda > 1e-4);
  * 低功率失稳窗口的上端 b_c: 取窗口内最后一个不稳定点 b1 与下一个稳定点 b2 的区间 [b1, b2],
    并用 Hamilton-Hopf 型阈值 growth^2 ∝ (b_c - b) 由最后两个不稳定点外推得到估计值 b_c_est (截断在 [b1, b2] 内);
  * 最大增长率及其 Im(lambda), 以及 b=7 处的增长率 (高功率端)。
若给出 --coarse (例如 hq_coarse), 同时比较两套网格的阈值区间。

输出 (results/analysis/<tag>/):
  stability_windows.csv      上述各量
  stability_map.png          (a) m=2 增长率 vs b; (b) m=1 增长率 vs b; (c) m>=3 高功率增长率 vs m

运行:  python a2_stability_map.py [--tag hq]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import style  # noqa: E402

TOL = 1e-4


def load(tag, kind, m):
    f = os.path.join(io_utils.PROJECT_ROOT, "results", "E3_stability", tag, f"{kind}_m{m}_spectrum.csv")
    if not os.path.exists(f):
        return None
    rows = io_utils.read_csv(f)
    b = np.array([float(r["b"]) for r in rows])
    g = np.array([float(r["max_growth"]) for r in rows])
    im = np.array([float(r["lam_im"]) for r in rows])
    o = np.argsort(b)
    return b[o], g[o], im[o]


def windows(b, g):
    """不稳定区间列表 [(b_start, b_end), ...] (按采样点)。"""
    out, cur = [], None
    for bi, gi in zip(b, g):
        if gi > TOL:
            cur = [bi, bi] if cur is None else [cur[0], bi]
        elif cur is not None:
            out.append(tuple(cur))
            cur = None
    if cur is not None:
        out.append(tuple(cur))
    return out


def low_window_threshold(b, g):
    """若最低的 b 点不稳定且之后出现稳定点: 返回 (b1, b2, b_c_est); 否则 None。"""
    if g[0] <= TOL:
        return None
    k = np.argmax(g <= TOL)
    if g[k] > TOL:                   # 整个范围都不稳定
        return None
    b1, b2 = b[k - 1], b[k]
    est = b2
    if k >= 2 and g[k - 2] > g[k - 1]:
        g0, g1 = g[k - 2] ** 2, g[k - 1] ** 2
        est = b1 + (b1 - b[k - 2]) * g1 / (g0 - g1)
        est = float(np.clip(est, b1, b2))
    return float(b1), float(b2), float(est)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    p.add_argument("--coarse", default="hq_coarse", help="粗网格对照 tag (不存在则跳过)")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(io_utils.PROJECT_ROOT, "results", "analysis", a.tag))
    rows = []
    data = {}
    for kind in style.SAMPLE_ORDER:
        for m in range(1, 7):
            d = load(a.tag, kind, m)
            if d is None:
                continue
            b, g, im = d
            data[(kind, m)] = d
            j = int(np.argmax(g))
            thr = low_window_threshold(b, g)
            dc = load(a.coarse, kind, m)
            thr_c = low_window_threshold(dc[0], dc[1]) if dc is not None else None
            g7 = float(g[np.argmin(abs(b - 7.0))]) if abs(b - 7.0).min() < 0.15 else float("nan")
            rows.append(dict(sample=kind, m=m, b_min=b[0], b_max=b[-1], n_spectra=len(b),
                             unstable_ranges=";".join(f"{x:.2f}-{y:.2f}" for x, y in windows(b, g)),
                             low_window_b1=thr[0] if thr else "", low_window_b2=thr[1] if thr else "",
                             b_c_est=thr[2] if thr else "", max_growth=g[j], b_at_max=b[j], im_at_max=im[j],
                             growth_b7=g7,
                             coarse_b1=thr_c[0] if thr_c else "", coarse_b2=thr_c[1] if thr_c else ""))
            if thr:
                print(f"{kind:5s} m={m}: low-power window up to b_c in [{thr[0]:.2f}, {thr[1]:.2f}] "
                      f"(est {thr[2]:.2f})" + (f"; coarse grid [{thr_c[0]:.2f}, {thr_c[1]:.2f}]" if thr_c else "")
                      + f"; max growth {g[j]:.4f} at b={b[j]:.2f} (Im {im[j]:.3f})")
    io_utils.write_csv(os.path.join(out, "stability_windows.csv"), rows)

    style.apply()
    fig, axs = style.plt.subplots(1, 3, figsize=(12.5, 3.7), gridspec_kw=dict(width_ratios=[1.25, 1.0, 0.9]))
    for ax, m, xmax in ((axs[0], 2, 3.4), (axs[1], 1, 1.6)):
        for kind in style.SAMPLE_ORDER:
            if (kind, m) not in data:
                continue
            b, g, im = data[(kind, m)]
            sel = b <= xmax
            st = style.sample_style(kind)
            ax.plot(b[sel], np.maximum(g[sel], 0.0), "-", color=st["color"], marker=st["marker"],
                    label=style.SAMPLE_LABEL[kind])
        ax.axhline(TOL, color=style.INK_SECONDARY, lw=0.8)
        ax.set_xlabel("propagation constant b")
        ax.set_ylabel("growth rate  max Re(lambda)")
        ax.set_title(f"({'a' if m == 2 else 'b'}) charge m={m}: low-power instability window", loc="left")
        ax.set_xlim(0.1, xmax)
        ax.set_ylim(bottom=0)
    axs[0].legend(loc="upper right")
    axs[0].text(0.62, 0.002, "pentagon = disk", fontsize=7.5, color=style.INK_SECONDARY)
    axs[1].legend(loc="upper right")
    ax = axs[2]
    for kind in style.SAMPLE_ORDER:
        ms = [m for m in range(3, 7) if (kind, m) in data]
        g7 = [float(data[(kind, m)][1][np.argmin(abs(data[(kind, m)][0] - 7.0))]) for m in ms]
        st = style.sample_style(kind)
        ax.plot(ms, g7, "-", color=st["color"], marker=st["marker"], label=style.SAMPLE_LABEL[kind])
    ax.text(4.5, 0.245, "n = 3-6 and disk coincide\n(beam far from the boundary)", ha="center", va="top",
            fontsize=7.5, color=style.INK_SECONDARY)
    ax.set_xticks([3, 4, 5, 6])
    ax.set_xlabel("topological charge m")
    ax.set_ylabel("growth rate at b = 7")
    ax.set_title("(c) m >= 3: unstable at all b", loc="left")
    ax.legend(loc="upper left")
    style.save(fig, os.path.join(out, "stability_map.png"))
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
