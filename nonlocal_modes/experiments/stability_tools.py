"""
stability_tools.py —— E3b、E5c 共用的稳定性工具
================================================
* low_power_shifts : 低功率失稳扫描用的 shift 集合 (omega = 0..1.2 密, 1.6..3.0 稀)
* spectrum_point   : 一个解的增长率、最不稳定本征值、模式的角向结构 (L、kappa、主导谐波)
* solve_between    : 两个已知解之间某个 b 处求解 (阈值二分用)
* bisect_thresholds: 稳定/不稳定交界处二分加密
* windows_and_thresholds: 由谱表得到失稳窗口与阈值区间
* plot_mode        : 模式图 (|u|、arg u、|v|、扰动后的光强)
"""

import json

import numpy as np

from common import pick

from nlmodes import diagnostics as dg
from nlmodes import newton, stability

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TOL = 1e-4


def low_power_shifts(args):
    om = np.arange(0.0, pick(args, quick=1.2, main=1.2, hq=1.2) + 1e-9, pick(args, quick=0.1, main=0.05, hq=0.05))
    tail = [1.6, 2.0, 2.5, 3.0]
    return [0.02 + 1j * w for w in list(om) + tail]


def spectrum_point(M, s, m, nsym, shifts, k):
    """一个解的稳定性: 增长率、最不稳定本征值及其模式的角向结构。"""
    st = stability.stability_arnoldi(M, s, shifts=shifts, k=k, return_mode=True)
    lam = st["eigenvalues"]
    j = stability.most_unstable_index(lam)            # λ 与 λ* 中取 Im λ >= 0 (L 的符号唯一)
    row = dict(b=s["b"], U=s["power"], growth=st["max_growth"], lam_re=float(lam[j].real),
               lam_im=float(abs(lam[j].imag)), unstable=st["max_growth"] > TOL,
               n_failed_shifts=st.get("n_failed_shifts", 0), time=round(st["time"], 1))
    mode = st.get("mode")
    if row["unstable"] and mode is not None:
        ac = stability.mode_angular_content(M, mode, m, nsym=nsym, w=s["w"])
        row.update(L=ac["L"], L_from_v=ac["L_from_v"], kappa=ac.get("kappa", ""), l_u=ac["l_u"], l_v=ac["l_v"],
                   top_u=json.dumps(ac["top_u"]), top_v=json.dumps(ac["top_v"]),
                   v_over_u=round(ac["norm_ratio"], 3))
    return row, mode


def solve_between(M, s1, s2, b):
    """在两个已知解之间的 b 处求解 (线性插值作初值)。"""
    t = (b - s1["b"]) / (s2["b"] - s1["b"])
    w2 = newton.align_phase(s2["w"], s1["w"])
    sol = newton.solve_robust(M, (1 - t) * s1["w"] + t * w2, b)
    return sol if sol["converged"] else None


def plot_mode(M, s, mode, row, fname, title):
    q = dg.field_grid(M, s["w"])
    u = dg.field_grid(M, mode["u"])
    v = dg.field_grid(M, mode["v"])
    eps = 0.25 * np.abs(q).max()
    pert = np.abs(q + eps * (u + np.conj(v))) ** 2
    ring = dg.ring_radius(q, M.grid, (0, 0))
    lim = min(2.2 * ring + 1.0, 0.5 * (M.grid.x[-1] - M.grid.x[0]))
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.2))
    panels = ((np.abs(u), "|u|", "magma"), (np.angle(u), "arg u", "twilight"), (np.abs(v), "|v|", "magma"),
              (pert, "|w + eps(u + v*)|^2", "inferno"))
    for ax, (f, t, cmap) in zip(axs, panels):
        ax.imshow(f, origin="lower", extent=M.grid.extent, cmap=cmap)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(t, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"{title}: b={s['b']:.3f}, lambda={row['lam_re']:.4f}+{row['lam_im']:.4f}i, "
                 f"L={row.get('L')} kappa={row.get('kappa')}  u{row.get('top_u')} v{row.get('top_v')}", fontsize=8)
    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    plt.close(fig)


def bisect_thresholds(M, rows, pts, modes, m, nsym, shifts, k, n_bisect, verbose=True):
    """在相邻两个谱点稳定性不同处二分 n_bisect 次; 新点追加到 rows, 模式存入 modes。"""
    coarse = sorted(zip([r["b"] for r in rows], [r["unstable"] for r in rows], pts), key=lambda t: t[0])
    for (b1, u1, s1), (b2, u2, s2) in zip(coarse[:-1], coarse[1:]):
        if u1 == u2:
            continue
        lo, hi = (s1, u1), (s2, u2)
        for _ in range(n_bisect):
            bm = 0.5 * (lo[0]["b"] + hi[0]["b"])
            sm = solve_between(M, lo[0], hi[0], bm)
            if sm is None:
                break
            row, mode = spectrum_point(M, sm, m, nsym, shifts, k)
            row["refined"] = True
            rows.append(row)
            modes[round(sm["b"], 6)] = (sm, mode)
            if verbose:
                print(f"   bisect b={bm:.4f} growth={row['growth']:.3e} L={row.get('L', '')}", flush=True)
            if row["unstable"] == lo[1]:
                lo = (sm, row["unstable"])
            else:
                hi = (sm, row["unstable"])
    rows.sort(key=lambda r: r["b"])


def windows_and_thresholds(rows):
    thresholds, windows, cur = [], [], None
    for r1, r2 in zip(rows[:-1], rows[1:]):
        if r1["unstable"] != r2["unstable"]:
            thresholds.append((r1["b"], r2["b"], "stable->unstable" if r2["unstable"] else "unstable->stable"))
    for r in rows:
        if r["unstable"]:
            cur = cur or []
            cur.append(r)
        elif cur:
            windows.append(cur)
            cur = None
    if cur:
        windows.append(cur)
    return windows, thresholds
