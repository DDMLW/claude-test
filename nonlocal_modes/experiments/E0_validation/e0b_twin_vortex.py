"""
E0b —— 复现导师论文: Ye, Kartashov, Hu, Torner, "Twin-vortex solitons in nonlocal nonlinear media",
        Opt. Lett. 35, 628 (2010)
=================================================================================================
目的: 用本代码复现已发表结果, 证明代码可靠 (审稿人一定会问)。坐标对应: 文中 eta -> x, zeta -> y,
样品 d_eta x d_zeta = 20 x 40 (x 方向是短轴)。

A) 文中 Fig.1: (m, b) = (1,2) (2,2) (3,4.64) (4,5) (5,5) (6,6) 的振幅/相位图
B) 文中 Fig.2(a), 3(b): m=2 解族 U(b) 与奇点间距 D_s(b)
C) 文中 Fig.2(b), 3(a): b=2, m=2 时积分宽度 W_eta, W_zeta 与 D_s 随 d_eta 的变化 (d_zeta=40)
D) 文中 Fig.4: 加噪声传播 (m=1,b=2), (m=2,b=2), (m=3,b=6); 以及 20x40 的精确解放进 21x39 样品中传播
E) 线性稳定性 (文中只做了传播, 这里补充本征值, 是对原文的加强)

运行:
  python e0b_twin_vortex.py --quick     (约 5-10 分钟, 粗网格 h=0.4, 短传播)
  python e0b_twin_vortex.py             (约 40-90 分钟, h=0.2, 传播到 xi=300)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (finalize_args, param_continuation, pick, pick_closest, rect_state_via_square,  # noqa: E402
                    robust_family, save_family, standard_args, transfer_solution)

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, propagate, stability  # noqa: E402
from nlmodes.geometry import rectangle_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FIG1_CASES = [(1, 2.0), (2, 2.0), (3, 4.64), (4, 5.0), (5, 5.0), (6, 6.0)]


def lg(model, m):
    return lambda W: guesses.lg_vortex(model, m, W)


def stability_summary(model, sol, shifts, k=12):
    st = stability.stability_arnoldi(model, sol, shifts=shifts, k=k)
    lam = st["eigenvalues"]
    k = int(np.argmax(lam.real))
    return dict(max_growth=st["max_growth"], lam_re=float(lam[k].real), lam_im=float(lam[k].imag))


def main():
    p = standard_args(__doc__)
    p.add_argument("--zmax", type=float, default=None, help="传播距离 (默认 300, quick 30, hq 1000)")
    args = finalize_args(p, "E0_validation")
    zmax = args.zmax or pick(args, quick=30.0, main=300.0, hq=1000.0)
    h = args.h_large                      # 20x40 矩形是“大样品”, hq 档用 h_large
    db = 0.25 * args.db_factor
    out, dat = args.res_dir, args.data_dir
    M = NonlocalModel(rectangle_sample(20, 40, h))
    print(M.describe())
    # hq: 更密、更宽的 shift 扫描 (omega 到 2, 每个 shift 求 16 个本征值), 避免漏掉虚部较大的不稳定模
    shifts = pick(args, quick=(0.01, 0.01 + 0.3j, 0.01 + 0.6j, 0.01 + 1.0j),
                  main=(0.01, 0.01 + 0.3j, 0.01 + 0.6j, 0.01 + 1.0j),
                  hq=tuple(0.01 + 1j * om for om in np.arange(0.0, 2.0 + 1e-9, 0.25)))
    kst = pick(args, quick=12, main=12, hq=16)

    # ------------------------------------------------------------------ A + E
    print("\n[A] Fig.1 cases + linear stability")
    fig1_rows = []
    fig1_sols = {}
    for m, bt in FIG1_CASES:
        if args.quick and m > 4:
            continue
        sols, info = robust_family(M, lg(M, m), bt, bt + 1.5, db=db, b_start=bt + 1.0, verbose=False)
        print(f"  m={m}: b in [{info.get('b_min')}, {info.get('b_max')}]  ({info['stop_down']})")
        s = pick_closest(sols, "b", bt)
        if s is None or abs(s["b"] - bt) > 1e-6:
            print(f"   direct solve failed for m={m}; trying aspect-ratio continuation from the 40x40 square")
            Mr, s = rect_state_via_square(h, m, bt)
            if s is None:
                print(f"   !! could not reach b={bt} for m={m}")
                continue
            s["w"] = transfer_solution(Mr, s["w"], M)
            s = newton.solve_robust(M, s["w"], bt)
            if not s["converged"]:
                print(f"   !! could not reach b={bt} for m={m}")
                continue
        fig1_sols[m] = s
        plotting.plot_mode(M, s, os.path.join(out, f"e0b_fig1_m{m}_b{bt:g}.png"),
                           title=f"20x40 sample, m={m}, b={bt:g}, U={s['power']:.3f}")
        io_utils.save_solution(os.path.join(dat, f"e0b_m{m}_b{bt:g}.npz"), M, s, extra=dict(m=m))
        row = dg.mode_summary(M, s, nsym=2, m=m)
        row.update(m=m, **stability_summary(M, s, shifts, kst))
        fig1_rows.append(row)
        print(f"   b={bt}: U={s['power']:.4f} central={row['central_charge']} sats={row['satellites']}"
              f"  max Re(lambda)={row['max_growth']:.2e}")
    io_utils.write_csv(os.path.join(out, "e0b_fig1_table.csv"), fig1_rows)

    # ------------------------------------------------------------------ B
    print("\n[B] m=2 family in 20x40: U(b), D_s(b)")
    b_lo = 0.5 if args.quick else 0.1
    sols, info = robust_family(M, lg(M, 2), b_lo, 6.0, db=db, b_start=4.0, verbose=False)
    print(f"  b range reached: [{info['b_min']:.3f}, {info['b_max']:.3f}];  lower end: {info['stop_down']}")
    rows = save_family(M, sols, out, dat, "e0b_m2_20x40", nsym=2, m=2, plot_bs=[1.0, 2.0, 4.0])
    for r in rows:
        # D_s: 两个卫星都被网格检测到时用检测值, 否则 (卫星落在中心 1.25h 以内) 用谐波拟合值
        r["D_s_detected"] = 2 * r["satellite_r_mean"]
        r["D_s_fit"] = 2 * r["r_s_fit"]
        r["D_s"] = 2 * (r["satellite_r_mean"] if r["n_satellites"] == 2 else r["r_s_fit"])
    io_utils.write_csv(os.path.join(out, "e0b_m2_20x40_summary.csv"), rows)
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.4))
    axs[0].plot([r["b"] for r in rows], [r["U"] for r in rows], "o-", ms=3)
    axs[0].set_xlabel("b")
    axs[0].set_ylabel("U")
    axs[0].set_title("Fig.2(a): m=2, 20x40")
    axs[1].plot([r["b"] for r in rows], [r["D_s"] for r in rows], "o-", ms=3)
    axs[1].set_xlabel("b")
    axs[1].set_ylabel("D_s")
    axs[1].set_title("Fig.3(b): d_eta=20")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "e0b_fig2a_fig3b.png"), dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------------ C
    print("\n[C] vary d_eta at b=2 (d_zeta=40)")
    s20 = fig1_sols.get(2)
    up = [25, 30, 35, 40] if not args.quick else [30, 40]
    down = [17.5, 15, 12.5, 10] if not args.quick else [15]
    deta_rows = []

    def build(d):
        return NonlocalModel(rectangle_sample(d, 40, h))

    def solve_b2(model, w0, prev):
        return newton.solve_robust(model, w0, 2.0)

    scan = [(20, M, s20)] if s20 is not None else []
    if s20 is not None:
        scan += param_continuation(build, 20, M, s20, up, solve_b2)
        scan += param_continuation(build, 20, M, s20, down, solve_b2)
    for d, Md, sol in scan:
        row = dg.mode_summary(Md, sol, nsym=2, m=2)
        row.update(d_eta=d, D_s_detected=2 * row["satellite_r_mean"], D_s_fit=2 * row["r_s_fit"],
                   D_s=2 * (row["satellite_r_mean"] if row["n_satellites"] == 2 else row["r_s_fit"]))
        deta_rows.append(row)
        if d in (10, 20, 40):
            plotting.plot_mode(Md, sol, os.path.join(out, f"e0b_deta{d:g}_m2_b2.png"),
                               title=f"{d:g}x40 sample, m=2, b=2")
        print(f"   d_eta={d:5.1f}: U={sol['power']:.4f}  W_eta={row['Wx']:.4f}  W_zeta={row['Wy']:.4f}  "
              f"central={row['central_charge']}  D_s={row['D_s']:.4f}")
    deta_rows.sort(key=lambda r: r["d_eta"])
    io_utils.write_csv(os.path.join(out, "e0b_vs_deta.csv"), deta_rows)
    if deta_rows:
        fig, axs = plt.subplots(1, 2, figsize=(8, 3.4))
        dd = [r["d_eta"] for r in deta_rows]
        axs[0].plot(dd, [r["Wx"] for r in deta_rows], "o-", label="W_eta (x)")
        axs[0].plot(dd, [r["Wy"] for r in deta_rows], "s-", label="W_zeta (y)")
        axs[0].set_xlabel("d_eta")
        axs[0].legend()
        axs[0].set_title("Fig.2(b): m=2, b=2")
        axs[1].plot(dd, [r["D_s"] for r in deta_rows], "o-")
        axs[1].set_xlabel("d_eta")
        axs[1].set_ylabel("D_s")
        axs[1].set_title("Fig.3(a): b=2")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "e0b_fig2b_fig3a.png"), dpi=130)
        plt.close(fig)

    # ------------------------------------------------------------------ D
    print(f"\n[D] noisy propagation up to xi={zmax:g}")
    snap = [0, zmax / 4, zmax / 2, 3 * zmax / 4, zmax]
    prop_rows = []
    for m, bt in [(1, 2.0), (2, 2.0), (3, 6.0)]:
        s = fig1_sols.get(m)
        if s is None:
            sols_m, _ = robust_family(M, lg(M, m), bt, bt + 1.5, db=db, b_start=bt + 1.0, verbose=False)
            s = pick_closest(sols_m, "b", bt)
        if s is None:
            continue
        q0 = propagate.add_noise(s["w"], 0.01, seed=1)
        res = propagate.propagate(M, q0, zmax, dz=args.dz, diag_every=25, snapshot_zs=snap, w_ref=s["w"],
                                  nsym=2, verbose=True)
        tag = f"e0b_prop_m{m}_b{s['b']:g}"
        plotting.plot_snapshots(M, res["snapshots"], os.path.join(out, tag + "_snapshots.png"),
                                zoom=7, title=f"m={m}, b={s['b']:g}, noise var 0.01", show_phase=True)
        plotting.plot_timeseries(res, os.path.join(out, tag + "_timeseries.png"), title=tag)
        prop_rows.append(dict(m=m, b=s["b"], zmax=zmax, final_deviation=float(res["deviation"][-1]),
                              max_deviation=float(res["deviation"].max()),
                              power_drift=float(np.ptp(res["power"]) / res["power"][0])))
        print(f"   m={m}: final deviation {res['deviation'][-1]:.3e}  (max {res['deviation'].max():.3e})")

    # 形状扰动: 20x40 的精确解放进 21x39 样品
    s2 = fig1_sols.get(2)
    if s2 is not None:
        Mp = NonlocalModel(rectangle_sample(21, 39, h))
        saved = dict(x=M.grid.x, y=M.grid.y, h=M.grid.h, q=dg.field_grid(M, s2["w"]))
        q0 = io_utils.solution_on_model(saved, Mp)
        res = propagate.propagate(Mp, q0, zmax, dz=args.dz, diag_every=25, snapshot_zs=snap,
                                  track_singularities=True, verbose=True)
        angles = []
        for z, sg in zip(res["z"], res["singularities"]):
            pos = [(x, y) for (x, y, qq) in sg if qq > 0 and np.hypot(x, y) < 4]
            if len(pos) == 2:
                (x1, y1), (x2, y2) = pos
                ang = np.degrees(np.arctan2(y2 - y1, x2 - x1)) % 180
                angles.append(dict(z=z, angle_deg=ang, separation=float(np.hypot(x2 - x1, y2 - y1))))
        io_utils.write_csv(os.path.join(out, "e0b_shape_perturbation_angle.csv"), angles)
        plotting.plot_snapshots(Mp, res["snapshots"], os.path.join(out, "e0b_shape_perturbation_snapshots.png"),
                                zoom=7, title="20x40 twin vortex launched into 21x39", show_phase=True)
        if angles:
            fig, ax = plt.subplots(figsize=(5, 3.2))
            ax.plot([a["z"] for a in angles], [a["angle_deg"] for a in angles])
            ax.set_xlabel("xi")
            ax.set_ylabel("angle of singularity line (deg, mod 180)")
            fig.tight_layout()
            fig.savefig(os.path.join(out, "e0b_shape_perturbation_angle.png"), dpi=130)
            plt.close(fig)
        prop_rows.append(dict(m=2, b=2.0, zmax=zmax, case="20x40 solution in 21x39 sample",
                              power_drift=float(np.ptp(res["power"]) / res["power"][0])))
    io_utils.write_csv(os.path.join(out, "e0b_propagation.csv"), prop_rows)
    print("\nE0b done.")


if __name__ == "__main__":
    main()
