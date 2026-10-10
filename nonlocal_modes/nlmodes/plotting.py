"""
plotting.py —— 标准作图 (所有图都保存为 PNG, 不弹窗)
=====================================================
图中文字统一用英文 (论文用图), 避免中文字体问题。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import diagnostics as dg  # noqa: E402


def _boundary(ax, model, color="w", lw=0.8):
    g = model.grid
    ax.contour(g.X, g.Y, model.opt_mask.astype(float), levels=[0.5], colors=color, linewidths=lw)


def _zoom_box(model, w, zoom):
    if zoom is None:
        g = model.grid
        return g.x[0], g.x[-1], g.y[0], g.y[-1]
    xc, yc = dg.center_of_mass(model, w)
    if zoom == "auto":
        Wx, Wy = dg.rms_widths(model, w)
        zoom = 3.5 * max(Wx, Wy, 1.0)
    return xc - zoom, xc + zoom, yc - zoom, yc + zoom


def mark_singularities(ax, sings, size=30):
    for (x, y, q) in sings:
        ax.scatter([x], [y], s=size, facecolors="none",
                   edgecolors="red" if q > 0 else "cyan", linewidths=1.5)


def singularities_for_plot(q, model, center, rel_threshold=1e-4):
    """
    画图用: 中心 1.25h 以外的单位奇点 (方格法) + 中心荷 (圆周绕数, 避免多重荷的混叠)。
    返回 (sings, central_charge)。
    """
    g = model.grid
    rc = 1.25 * g.h
    reg = dg.beam_region(q, g, rel_threshold=rel_threshold) & model.opt_mask
    sings = dg.find_singularities(q, g, reg, exclude_center=center, exclude_radius=rc)
    central = int(round(dg.winding_on_circle(q, g, center, 0.5 * g.h)))
    return sings, central


def annotate_center(ax, center, central, fontsize=8):
    if central != 0:
        ax.plot(center[0], center[1], "x", color="red" if central > 0 else "cyan", ms=6, mew=1.5)
        ax.annotate(f"{central:+d}", center, color="w", fontsize=fontsize, xytext=(4, 4),
                    textcoords="offset points")


def plot_mode(model, sol, fname, title=None, zoom="auto", show_index=True, show_sings=True,
              center=(0.0, 0.0)):
    """
    三联图: |w|, arg(w)(红圈 +1、青圈 -1 为单位奇点; 中心的叉号与数字为圆周绕数给出的中心荷),
    折射率 n (整个样品)。center=None 时以光束质心为中心 (表面态等不在样品中心的态)。
    """
    w = sol["w"]
    q = dg.field_grid(model, w)
    g = model.grid
    if center is None:
        center = dg.center_of_mass(model, w)
    ncol = 3 if show_index else 2
    fig, axs = plt.subplots(1, ncol, figsize=(4.2 * ncol, 4))
    x0, x1, y0, y1 = _zoom_box(model, w, zoom)
    sings, central = [], 0
    if show_sings:
        sings, central = singularities_for_plot(q, model, center)
    im = axs[0].imshow(np.abs(q), origin="lower", extent=g.extent, cmap="inferno")
    axs[0].set_title("|q|")
    fig.colorbar(im, ax=axs[0], fraction=0.046)
    ph = np.where(np.abs(q) > 1e-6 * np.abs(q).max(), np.angle(q), np.nan)
    im = axs[1].imshow(ph, origin="lower", extent=g.extent, cmap="twilight", vmin=-np.pi, vmax=np.pi)
    axs[1].set_title("arg q  (red: +1, cyan: -1)")
    fig.colorbar(im, ax=axs[1], fraction=0.046)
    mark_singularities(axs[1], sings)
    if show_sings:
        annotate_center(axs[1], center, central)
    for ax in axs[:2]:
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        _boundary(ax, model)
    if show_index:
        n_t = sol.get("n_t")
        if n_t is None:
            n_t = model.thermal_index(np.abs(w) ** 2)
        im = axs[2].imshow(model.index_grid(n_t), origin="lower", extent=g.extent, cmap="viridis")
        axs[2].set_title("n (whole sample)")
        fig.colorbar(im, ax=axs[2], fraction=0.046)
        _boundary(axs[2], model, color="k")
    for ax in axs:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.suptitle(title or f"b={sol['b']:.4g}, U={sol['power']:.4g}")
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)
    return sings


def plot_family(rows, fname, x="b", ys=("U",), title=None, logy=False):
    fig, axs = plt.subplots(1, len(ys), figsize=(4.2 * len(ys), 3.6))
    axs = np.atleast_1d(axs)
    for ax, y in zip(axs, ys):
        xs = [float(r[x]) for r in rows]
        vs = [float(r[y]) for r in rows]
        ax.plot(xs, vs, "o-", ms=3)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if logy:
            ax.set_yscale("log")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def plot_spectrum(lams, fname, title=None):
    fig, ax = plt.subplots(figsize=(4.5, 4))
    ax.plot(lams.real, lams.imag, ".", ms=3)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("Re lambda")
    ax.set_ylabel("Im lambda")
    ax.set_title(title or f"max Re = {np.max(lams.real):.2e}")
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def plot_snapshots(model, snapshots, fname, zoom=None, center=(0.0, 0.0), title=None, show_phase=False):
    """传播快照: 第一行 |q|, (可选)第二行相位并标出奇点。"""
    zs = sorted(snapshots)
    g = model.grid
    nrow = 2 if show_phase else 1
    fig, axs = plt.subplots(nrow, len(zs), figsize=(2.6 * len(zs), 2.6 * nrow), squeeze=False)
    for k, z in enumerate(zs):
        q = dg.field_grid(model, snapshots[z])
        ax = axs[0, k]
        ax.imshow(np.abs(q), origin="lower", extent=g.extent, cmap="inferno")
        ax.set_title(f"xi={z:g}", fontsize=9)
        if zoom:
            ax.set_xlim(center[0] - zoom, center[0] + zoom)
            ax.set_ylim(center[1] - zoom, center[1] + zoom)
        ax.set_xticks([])
        ax.set_yticks([])
        if show_phase:
            ax = axs[1, k]
            ph = np.where(np.abs(q) > 1e-3 * np.abs(q).max(), np.angle(q), np.nan)
            ax.imshow(ph, origin="lower", extent=g.extent, cmap="twilight", vmin=-np.pi, vmax=np.pi)
            sg, central = singularities_for_plot(q, model, center, rel_threshold=1e-3)
            mark_singularities(ax, sg, size=15)
            annotate_center(ax, center, central, fontsize=7)
            if zoom:
                ax.set_xlim(center[0] - zoom, center[0] + zoom)
                ax.set_ylim(center[1] - zoom, center[1] + zoom)
            ax.set_xticks([])
            ax.set_yticks([])
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    plt.close(fig)


def plot_timeseries(res, fname, keys=("power", "peak", "deviation"), title=None):
    keys = [k for k in keys if k in res and len(res[k])]
    fig, axs = plt.subplots(1, len(keys), figsize=(4 * len(keys), 3.2))
    axs = np.atleast_1d(axs)
    for ax, k in zip(axs, keys):
        ax.plot(res["z"], res[k])
        ax.set_xlabel("xi")
        ax.set_ylabel(k)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    plt.close(fig)
