"""
fieldplot.py —— 论文场图的公用函数 (a9 的 Fig.2 图谱、a7 的 Fig.6 双涡旋场图)
==========================================================================
读取实验脚本导出的小 npz (E2c: results/E2_transmutation_atlas/<tag>/fields/; E5d --fields_only:
results/E5_thermal_reconfiguration/<tag>/fields/), 键: x, y (格心坐标), psi (复光场), h, 以及若干标量元数据。
着色: 色相 = 相位 (HSV), 亮度 = sqrt(|psi| / max|psi|) —— 不能用 twilight 之类本身明暗变化的循环色图,
否则相位与振幅都会表现为明暗。奇点: 方格法 (中心 1.25h 以外, 只在光束区内), 中心荷用 0.5h 圆周绕数 (numerics.md §6)。
"""

import numpy as np

from nlmodes import diagnostics as dg
from nlmodes.geometry import Grid

import style


def load_field(f):
    """返回 (Grid, psi, meta); meta 是 npz 中全部标量 (b、U、h、Theta、m ...)。"""
    d = np.load(f, allow_pickle=False)
    g = Grid.__new__(Grid)
    g.x, g.y, g.h = d["x"].astype(float), d["y"].astype(float), float(d["h"])
    g.X, g.Y = np.meshgrid(g.x, g.y)
    g.shape = g.X.shape
    g.R, g.TH = np.hypot(g.X, g.Y), np.arctan2(g.Y, g.X)
    meta = {k: d[k].item() for k in d.files if d[k].ndim == 0}
    return g, d["psi"].astype(complex), meta


def domain_color(q):
    from matplotlib.colors import hsv_to_rgb
    hsv = np.stack([(np.angle(q) + np.pi) / (2 * np.pi), np.full(q.shape, 0.85),
                    np.sqrt(np.abs(q) / np.abs(q).max())], axis=-1)
    return hsv_to_rgb(hsv)


def draw_field(ax, g, q, half=None, ms=4.5):
    """在 ax 上画相位-振幅图并标出 +1 (白圈) / -1 (黑圈) 奇点; half 给定时只显示 |x|,|y| <= half。
    返回 (奇点列表 [(x, y, 荷)], 中心荷)。"""
    ext = [g.x[0] - g.h / 2, g.x[-1] + g.h / 2, g.y[0] - g.h / 2, g.y[-1] + g.h / 2]
    ax.imshow(domain_color(q), origin="lower", extent=ext, interpolation="bilinear")
    # 只在光束区内找奇点: 平滑光强 > 1e-2 倍峰值 (裁剪框可能越过多边形的边, 边外 q = 0 会给出成片的伪奇点)
    reg = dg.beam_region(q, g, rel_threshold=1e-2) & (np.abs(q) > 0)
    sings = dg.find_singularities(q, g, reg, exclude_center=(0.0, 0.0), exclude_radius=1.25 * g.h)
    for (x, y, c) in sings:
        if half is None or (abs(x) <= half and abs(y) <= half):
            ax.plot(x, y, "o", ms=ms, mfc="none", mec="white" if c > 0 else "black", mew=1.2)
    central = int(round(dg.winding_on_circle(q, g, (0.0, 0.0), 0.5 * g.h)))
    if half is not None:
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return sings, central


def phase_bar(ax):
    """横向的相位色条 (亮度 1), 刻度 -pi, 0, pi。"""
    ph = np.linspace(-np.pi, np.pi, 256)
    ax.imshow(domain_color(np.exp(1j * ph))[None, :, :], aspect="auto", extent=[-np.pi, np.pi, 0, 1])
    ax.set_yticks([])
    ax.set_xticks([-np.pi, 0, np.pi])
    ax.set_xticklabels(["-pi", "0", "pi"], fontsize=7)
    ax.set_xlabel("phase", fontsize=7, color=style.INK_SECONDARY)
    ax.grid(False)
