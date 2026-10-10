"""
E0c —— 复现导师论文: Ye, Kartashov, Torner, "Nonlocal surface dipoles and vortices", PRA 77, 033829 (2008)
======================================================================================================
几何: 热介质 -L/2 < x < 0, |y| < L/2 (L=40), 线性介质 x > 0 (n_d = -100), 界面 x=0 绝热。
(文中 eta -> x 是垂直于界面的坐标, zeta -> y 平行于界面。)

A) b=3 的表面偶极 (节线垂直于界面, 即两个光斑沿 y 排列)、表面涡旋、双涡旋束缚态 (文中 Fig.1-3)
B) 偶极 U(b) (文中 Fig.4a) 与双涡旋束缚态 U(b) (Fig.4b)
C) 涡旋的非正则相位: 绕奇点一周的相位 theta(phi) (文中 Fig.4c)
D) 加噪声 (方差 0.01) 传播 (文中 Fig.5) + 线性稳定性 (原文没有, 作为加强)

运行:
  python e0c_surface_solitons.py --quick   (约 5-10 分钟)
  python e0c_surface_solitons.py           (约 1-2 小时)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, initial_solution, save_family, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, propagate, stability  # noqa: E402
from nlmodes.geometry import surface_interface_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _interface_factor(M, W):
    """线性介质 (n_d=-100) 强烈排斥光, 表面态在界面 x=0 处几乎为零: 初值乘以 max(-x,0)/W。"""
    return np.clip(-M.xo, 0, None) / W


def dipole_shape(M, x0=-1.5, d=1.6):
    """节线垂直于界面的偶极: 两个反相光斑位于 (x0, +-d)。"""
    def f(W):
        return (_interface_factor(M, W) * (np.exp(-((M.xo - x0) ** 2 + (M.yo - d) ** 2) / W ** 2)
                                           - np.exp(-((M.xo - x0) ** 2 + (M.yo + d) ** 2) / W ** 2))
                ).astype(complex)
    return f


def vortex_shape(M, x0=-2.0):
    def f(W):
        z = (M.xo - x0) + 1j * M.yo
        return _interface_factor(M, W) * (z / W) * np.exp(-np.abs(z) ** 2 / W ** 2)
    return f


def bound_state_shape(M, x0=-3.0, d=1.6):
    """两个同号单位涡旋沿 x (垂直界面) 排列, 嵌在一个沿 x 拉长的包络里。"""
    def f(W):
        z1 = (M.xo - (x0 + d)) + 1j * M.yo
        z2 = (M.xo - (x0 - d)) + 1j * M.yo
        env = np.exp(-(((M.xo - x0) / (1.6 * W)) ** 2 + (M.yo / W) ** 2))
        return _interface_factor(M, W) * z1 * z2 / W ** 2 * env
    return f


def search_complex_state(M, U, pairs, kappas, accept, x0=-2.5, w0=2.0):
    """
    自洽线性模式初值 + 定功率牛顿: 对模式对 (i, j) 与混合比 kappa 扫描, 场 = psi_i + i*kappa*psi_j。
    accept(sol, row) 判断收敛解是否是想要的态 (例如: 两个同号奇点 = 涡旋束缚态)。
    表面涡旋是“非正则”的 (两个偶极模式权重不同), 所以必须扫描 kappa (见 docs/numerics.md)。
    """
    I0 = np.exp(-((M.xo - x0) ** 2 + M.yo ** 2) / w0 ** 2) * (M.xo < 0)
    for (i, j) in pairs:
        for kap in kappas:
            w, b = guesses.scf_guess(M, U, lambda v, e, i=i, j=j, kap=kap: v[:, i] + 1j * kap * v[:, j], I0)
            sol = newton.solve_robust(M, w, b, fixed_power=U, maxit=30)
            if sol["converged"]:
                row = dg.mode_summary(M, sol, center=None)
                if accept(sol, row):
                    print(f"   found with modes ({i},{j}), kappa={kap}: b={sol['b']:.4f}")
                    return sol
    return None


def _nodal_counts(M, v, xc, off=0.7):
    """
    沿水平线 y=+off (平行于 x, 即垂直界面方向) 与竖直线 x=xc+off 上的变号次数, 用于识别 HG 型线性模式。
    取偏离对称轴的线, 避免奇宇称模式在对称轴上恰为 0 造成误判。
    """
    g = M.grid
    vg = M.to_grid(v)
    row = vg[np.argmin(np.abs(g.y - off)), :]
    col = vg[:, np.argmin(np.abs(g.x - (xc + off)))]

    def changes(line):
        a = line[np.abs(line) > 0.05 * np.abs(line).max()]
        return int(np.sum(np.diff(np.sign(a)) != 0))
    return changes(row), changes(col)


def search_bound_state(M, U=8.0, x0=-3.5, w0=2.5):
    """
    双涡旋束缚态 ≈ (沿 x 有两个节点的模式) + i kappa (x、y 各一个节点的模式)
    —— 正是两个同号奇点排在 x 轴 (垂直界面) 上的结构。每次自洽迭代都按节点结构重新挑选模式。
    """
    I0 = np.exp(-(((M.xo - x0) / (1.6 * w0)) ** 2 + (M.yo / w0) ** 2)) * (M.xo < 0)

    def pick(vecs, kind, xc):
        for k in range(vecs.shape[1]):
            if _nodal_counts(M, vecs[:, k], xc) == kind:
                return vecs[:, k]
        return None

    for real_kind in ((2, 0), (2, 2), (0, 2)):
        for kap in (0.5, 1.0, 2.0):
            def combine(vecs, vals, kap=kap, real_kind=real_kind):
                I = np.abs(vecs[:, 0]) ** 2
                xc = float((M.xo * I).sum() / I.sum())
                a, b_ = pick(vecs, real_kind, xc), pick(vecs, (1, 1), xc)
                if a is None or b_ is None:
                    return vecs[:, 1] + 1j * kap * vecs[:, 2]
                return a + 1j * kap * b_
            w, b = guesses.scf_guess(M, U, combine, I0, k=12)
            sol = newton.solve_robust(M, w, b, fixed_power=U, maxit=30)
            if sol["converged"]:
                row = dg.mode_summary(M, sol, center=None)
                if abs(abs(row["total_charge"]) - 2) < 0.1:
                    print(f"   bound state found ({real_kind}, kappa={kap}): b={sol['b']:.4f}, "
                          f"singularities {row['satellites']}")
                    return sol
    return None


def main():
    p = standard_args(__doc__)
    p.add_argument("--zmax", type=float, default=None, help="传播距离 (默认 300, quick 30)")
    args = finalize_args(p, "E0_validation")
    zmax = args.zmax or (30.0 if args.quick else 300.0)
    out, dat = args.res_dir, args.data_dir
    M = NonlocalModel(surface_interface_sample(L=40.0, nd=-100.0, h=args.h))
    print(M.describe())
    b0 = 3.0
    shifts = (0.01, 0.01 + 0.3j, 0.01 + 0.6j, 0.01 + 1.0j)
    cases = dict(dipole=(dipole_shape(M), True), vortex=(vortex_shape(M), False),
                 bound_state=(bound_state_shape(M), False))
    sols_b3 = {}
    rows = []
    print("\n[A] b=3 states")
    for name, (shape, real) in cases.items():
        s = initial_solution(M, shape, [b0, b0 + 0.5, b0 - 0.5, b0 + 1.0], real=real, verbose=True)
        if s is None and name == "vortex":
            print("   LG-type guesses failed; trying self-consistent linear-mode guesses ...")
            s = search_complex_state(M, 6.0, [(2, 1), (1, 2)], [0.3, 0.5, 0.7, 1.4, 2.0],
                                     lambda sol, row: abs(abs(row["total_charge"]) - 1) < 0.1)
        if s is None and name == "bound_state":
            print("   product guesses failed; trying self-consistent linear-mode guesses (nodal classification) ...")
            s = search_bound_state(M)
        if s is None:
            print(f"   !! {name}: no convergence")
            continue
        if abs(s["b"] - b0) > 1e-9:
            sl, _ = newton.continuation(M, s["w"], s["b"], b0, np.sign(b0 - s["b"]) * 0.25, real=real,
                                        verbose=False)
            s = sl[-1]
        sols_b3[name] = s
        io_utils.save_solution(os.path.join(dat, f"e0c_{name}_b3.npz"), M, s)
        plotting.plot_mode(M, s, os.path.join(out, f"e0c_{name}_b3.png"), zoom=8,
                           title=f"surface {name}, b=3, U={s['power']:.3f}", center=None)
        row = dg.mode_summary(M, s, center=None)
        st = stability.stability_arnoldi(M, s, shifts=shifts, k=12)
        row.update(state=name, max_growth=st["max_growth"])
        rows.append(row)
        print(f"   {name}: U={s['power']:.4f} xc={row['xc']:.3f} sings={row['satellites']} "
              f"total charge={row['total_charge']}  max Re(lambda)={st['max_growth']:.2e}")
    io_utils.write_csv(os.path.join(out, "e0c_b3_states.csv"), rows)

    print("\n[B] U(b) for dipole and vortex bound state (+ stability along the branch)")
    b_lo = 0.5 if args.quick else 0.1
    b_hi = 5.0 if args.quick else 8.0
    fam_rows = {}
    for name in ("dipole", "bound_state"):
        s = sols_b3.get(name)
        if s is None:
            continue
        real = cases[name][1]
        up, _ = newton.continuation(M, s["w"], b0, b_hi, 0.25, real=real, verbose=False)
        dn, info = newton.continuation(M, s["w"], b0, b_lo, -0.25, real=real, verbose=False)
        fam = sorted(dn[1:] + up, key=lambda x: x["b"])
        print(f"   {name}: b in [{fam[0]['b']:.3f}, {fam[-1]['b']:.3f}]  lower end: {info['stop_reason']}")
        rr = save_family(M, fam, out, dat, f"e0c_{name}_family", center=None)
        every = 2 if args.quick else 4
        for r, sol in zip(rr[::every], fam[::every]):
            r["max_growth"] = stability.stability_arnoldi(M, sol, shifts=shifts, k=12)["max_growth"]
        io_utils.write_csv(os.path.join(out, f"e0c_{name}_family_summary.csv"), rr)
        fam_rows[name] = rr
    if fam_rows:
        fig, axs = plt.subplots(1, len(fam_rows), figsize=(4.2 * len(fam_rows), 3.4))
        axs = np.atleast_1d(axs)
        for ax, (name, rr) in zip(axs, fam_rows.items()):
            ax.plot([r["b"] for r in rr], [r["U"] for r in rr], "-", color="0.6")
            st = [r for r in rr if "max_growth" in r]
            ax.plot([r["b"] for r in st if r["max_growth"] < 1e-4], [r["U"] for r in st if r["max_growth"] < 1e-4],
                    "o", color="C0", label="stable")
            ax.plot([r["b"] for r in st if r["max_growth"] >= 1e-4], [r["U"] for r in st if r["max_growth"] >= 1e-4],
                    "x", color="C3", label="unstable")
            ax.set_xlabel("b")
            ax.set_ylabel("U")
            ax.set_title(name)
            ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out, "e0c_fig4ab.png"), dpi=130)
        plt.close(fig)

    print("\n[C] noncanonical phase of the surface vortex (Fig.4c)")
    s = sols_b3.get("vortex")
    if s is not None:
        q = dg.field_grid(M, s["w"])
        sg = [x for x in dg.find_singularities(q, M.grid, dg.beam_region(q, M.grid) & M.opt_mask) if x[2] == 1]
        if sg:
            xs, ys, _ = sg[0]
            th = np.linspace(0, 2 * np.pi, 361)
            r0 = 1.0
            vals = dg.interp(q, M.grid, xs + r0 * np.cos(th), ys + r0 * np.sin(th))
            ph = np.unwrap(np.angle(vals))
            ph -= ph[0]
            io_utils.write_csv(os.path.join(out, "e0c_vortex_phase_contour.csv"),
                               [dict(phi=float(a), theta=float(b_)) for a, b_ in zip(th, ph)])
            fig, ax = plt.subplots(figsize=(4.5, 3.4))
            ax.plot(th, ph, label="surface vortex")
            ax.plot(th, th, "k:", label="canonical")
            ax.set_xlabel("azimuth phi around singularity")
            ax.set_ylabel("phase theta")
            ax.legend()
            fig.tight_layout()
            fig.savefig(os.path.join(out, "e0c_fig4c_phase.png"), dpi=130)
            plt.close(fig)

    print(f"\n[D] noisy propagation (variance 0.01) up to xi={zmax:g}")
    snap = [0, zmax / 4, zmax / 2, 3 * zmax / 4, zmax]
    prow = []
    for name, s in sols_b3.items():
        xc, yc = dg.center_of_mass(M, s["w"])
        res = propagate.propagate(M, propagate.add_noise(s["w"], 0.01, seed=2), zmax, dz=0.02, diag_every=25,
                                  snapshot_zs=snap, w_ref=s["w"], verbose=True)
        plotting.plot_snapshots(M, res["snapshots"], os.path.join(out, f"e0c_prop_{name}.png"), zoom=8,
                                center=(xc, yc), title=f"surface {name}, b=3", show_phase=True)
        plotting.plot_timeseries(res, os.path.join(out, f"e0c_prop_{name}_timeseries.png"), title=name)
        prow.append(dict(state=name, zmax=zmax, final_deviation=float(res["deviation"][-1]),
                         max_deviation=float(res["deviation"].max())))
        print(f"   {name}: final deviation {res['deviation'][-1]:.3e}")
    io_utils.write_csv(os.path.join(out, "e0c_propagation.csv"), prow)
    print("\nE0c done.")


if __name__ == "__main__":
    main()
