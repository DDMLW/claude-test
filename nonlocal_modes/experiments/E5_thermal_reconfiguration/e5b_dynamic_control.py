"""
E5b —— 沿传播方向动态改变边界温度图样: 拓扑结构的“开关”、“切换”与“绕转”
========================================================================
设定: 圆形样品 (R=10) 中先求出 C_3 加热图样 (T_b = Theta cos 3theta) 下 m=2 的嬗变定态
(中心 -1 + 3 个 +1 卫星), 再让它在下列沿 xi 变化的热边界中传播 (分段加热片沿样品长度方向排列,
热传导仍按横向稳态处理, 要求图样沿 xi 的变化远慢于横向热扩散 —— 见 README 中的适用条件):

  (i)   关断:   xi > xi_s 时 Theta -> 0             问: 嬗变是否可逆 (卫星收回中心)?
  (ii)  切换:   xi > xi_s 时 C_3 图样 -> C_4 图样   问: m=2 在 C_4 中 mu = 2 = n/2 (简并情形), 结构如何演化?
  (iii) 旋转:   theta0(xi) = Omega * xi             问: 卫星环能否绝热地跟随旋转? 临界 Omega?
                (卫星的世界线在 (x, y, xi) 空间中形成螺旋/辫子)

输出: control_<case>_snapshots.png, control_<case>_tracks.csv/png (卫星方位角/半径随 xi), control_summary.csv

运行:
  python e5b_dynamic_control.py --quick    (约 5-10 分钟)
  python e5b_dynamic_control.py            (约 1-2 小时)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, initial_solution, param_continuation, pick, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, propagate  # noqa: E402
from nlmodes.geometry import disk_sample, harmonic_multipole  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def pattern_field(R, h, n, Theta, theta0):
    """返回热介质点上由边界温度 Theta cos(n(theta - theta0)) 产生的折射率 (调和场)。"""
    s = disk_sample(R, h)
    s.set_boundary_temperature(harmonic_multipole(n, Theta, R, theta0))
    return NonlocalModel(s).n_bc


def satellite_tracks(res, rmax):
    rows = []
    for z, sg in zip(res["z"], res["singularities"]):
        for (x, y, c) in sg:
            r = np.hypot(x, y)
            if r < rmax:
                rows.append(dict(z=z, x=x, y=y, r=r, angle=float(np.arctan2(y, x)), charge=c))
    return rows


def main():
    p = standard_args(__doc__)
    p.add_argument("--R", type=float, default=10.0)
    p.add_argument("--Theta", type=float, default=4.0,
                   help="C_3 图样的边界温度幅度 (默认 4: h=0.2 时卫星环 r~0.6 > 2h, 可分辨; 1 时卫星在亚网格)")
    p.add_argument("--U", type=float, default=15.0)
    p.add_argument("--zmax", type=float, default=None, help="默认 300 (quick 40, hq 600)")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    R, h, Th, U = args.R, args.h, args.Theta, args.U
    zmax = args.zmax or pick(args, quick=40.0, main=300.0, hq=600.0)
    zs = zmax / 4                                        # 切换位置
    m = 2

    # ---- 起点: C_3 加热下 m=2 的定态 (从 Theta=0 沿 Theta 延拓)
    M0 = NonlocalModel(disk_sample(R, h))
    print("computing the start state (Theta=0 -> Theta by parameter continuation) ...", flush=True)
    first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [3.0, 3.5, 2.5], verbose=True)
    s0 = newton.solve_robust(M0, first["w"], first["b"], fixed_power=U)

    def build(T):
        s = disk_sample(R, h)
        s.set_boundary_temperature(harmonic_multipole(3, T, R))
        return NonlocalModel(s)

    scan = param_continuation(build, 0.0, M0, s0, [Th / 16, Th / 8, Th / 4, Th / 2, Th],
                              lambda model, w0, prev: newton.solve_robust(model, w0, prev["b"], fixed_power=U))
    if not scan or scan[-1][0] != Th:
        raise SystemExit("未能得到 C_3 加热下的定态, 请减小 --Theta")
    _, M3, st3 = scan[-1]
    row = dg.mode_summary(M3, st3, nsym=3, m=m)
    print(f"start state: b={st3['b']:.4f}, central (small circle) {row['W_small']:+.0f}, "
          f"satellites {row['n_satellites']} at r={row['satellite_r_mean']:.3f}")
    plotting.plot_mode(M3, st3, os.path.join(out, "control_start_state.png"),
                       title=f"m=2, C3 heating Theta={Th:g}, U={U:g}")

    H3 = pattern_field(R, h, 3, Th, 0.0)
    H3s = pattern_field(R, h, 3, Th, np.pi / 6)          # 旋转 pi/6 (C_3 图样的 sin 分量)
    H4 = pattern_field(R, h, 4, Th, 0.0)
    ring = dg.ring_radius(dg.field_grid(M3, st3["w"]), M3.grid, (0, 0), rmax=8.0)

    def rotating(Om):
        # cos(3(theta - Om z)) = cos(3 Om z) cos 3theta + sin(3 Om z) sin 3theta, 而 H3s 对应 sin 3theta
        return lambda z: np.cos(3 * Om * z) * H3 + np.sin(3 * Om * z) * H3s

    cases = {
        "switch_off": lambda z: H3 if z < zs else 0.0 * H3,
        "switch_C3_to_C4": lambda z: H3 if z < zs else H4,
        "rotate_slow": rotating(0.005 if not args.quick else 0.02),
        "rotate_fast": rotating(0.05 if not args.quick else 0.2),
    }
    summary = []
    snap = list(np.linspace(0, zmax, 6))
    for name, f in cases.items():
        res = propagate.propagate(M0, st3["w"], zmax, dz=args.dz, diag_every=25, snapshot_zs=snap, nbc_func=f,
                                  track_singularities=True, verbose=False)
        plotting.plot_snapshots(M0, res["snapshots"], os.path.join(out, f"control_{name}_snapshots.png"),
                                zoom=2 * ring + 1, title=name, show_phase=True)
        tr = satellite_tracks(res, 1.5 * ring)
        io_utils.write_csv(os.path.join(out, f"control_{name}_tracks.csv"), tr)
        fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.3))
        for t in tr:
            col = "red" if t["charge"] > 0 else "blue"
            axs[0].plot(t["z"], t["r"], ".", ms=2, color=col)
            axs[1].plot(t["z"], t["angle"], ".", ms=2, color=col)
        axs[0].set_xlabel("xi")
        axs[0].set_ylabel("r of singularities")
        axs[1].set_xlabel("xi")
        axs[1].set_ylabel("angle of singularities")
        fig.suptitle(name)
        fig.tight_layout()
        fig.savefig(os.path.join(out, f"control_{name}_tracks.png"), dpi=130)
        plt.close(fig)
        qf = dg.field_grid(M0, res["final"])
        summary.append(dict(case=name, zmax=zmax, central_final=round(dg.winding_on_circle(qf, M0.grid, (0, 0),
                                                                                            0.5 * h), 2),
                            power_drift=float(np.ptp(res["power"]) / res["power"][0])))
        print(f"  {name}: final central charge {summary[-1]['central_final']:+.2f}  [{res['time']:.0f}s]")
    io_utils.write_csv(os.path.join(out, "control_summary.csv"), summary)
    print("\nE5b done.")


if __name__ == "__main__":
    main()
