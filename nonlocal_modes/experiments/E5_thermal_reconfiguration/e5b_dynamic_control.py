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

输出: control_<case>_snapshots.png, control_<case>_tracks.csv/png (卫星方位角/半径随 xi),
      control_<case>_timeseries.csv/png (C_n 非 mu 权重、与起点的偏离、中心荷), control_summary.csv
      (--n 不等于 3 时文件名前缀为 control_n<n>_T<Theta>_, 汇总为 control_n<n>_T<Theta>_summary.csv)

2026-10-12 推广 (--n, --m): 第二轮 E5c 发现 C_3 加热的 m=2 起点 (Theta=3, U=15) 线性不稳定 (增长率 0.105),
E5b 原有四条轨迹在 xi≈26-160 之后由不稳定性主导; 而 C_2 (四极) 加热下的 m=2 "双涡旋" 在云端 h=0.2 检查中稳定,
且卫星间距可调 (r_s ∝ Theta^(1/2), Theta=0.5/1 时 r_s=0.67/0.95)。--n 2 即用四极加热做同样的关断/旋转实验
(n=2 时没有 C_3 -> C_4 切换这一情形)。

运行:
  python e5b_dynamic_control.py --quick    (约 5-10 分钟)
  python e5b_dynamic_control.py            (约 1-2 小时)
  python e5b_dynamic_control.py --hq --n 2 --Theta 0.5 --job e5b_n2_T0.5     (四极加热双涡旋)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, initial_solution, param_continuation, pick, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, propagate, theory  # noqa: E402
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
    p.add_argument("--n", type=int, default=3, help="加热图样的阶数 (3: 原 C_3 情形; 2: 四极加热的双涡旋)")
    p.add_argument("--m", type=int, default=2, help="拓扑荷")
    p.add_argument("--Theta", type=float, default=3.0,
                   help="C_3 图样的边界温度幅度 (默认 3: 卫星环 r~0.44, h=0.2 时 >2h、h=0.1 时 ~4h, 可分辨; "
                        "1 时卫星在亚网格; U=15 时 Theta >~ 3.8 涡旋会与样品边缘的热点模式杂化, 不要超过 3.5)")
    p.add_argument("--U", type=float, default=15.0)
    p.add_argument("--zmax", type=float, default=None, help="默认 300 (quick 40, hq 600)")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    R, h, Th, U, n, m = args.R, args.h, args.Theta, args.U, args.n, args.m
    zmax = args.zmax or pick(args, quick=40.0, main=300.0, hq=600.0)
    zs = zmax / 4                                        # 切换位置
    pre = "control_" if n == 3 else f"control_n{n}_T{Th:g}_"     # n=3 保持原文件名

    # ---- 起点: C_n 加热下 m 荷涡旋的定态 (从 Theta=0 沿 Theta 延拓)
    M0 = NonlocalModel(disk_sample(R, h))
    print("computing the start state (Theta=0 -> Theta by parameter continuation) ...", flush=True)
    first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [3.0, 3.5, 2.5], verbose=True)
    s0 = newton.solve_robust(M0, first["w"], first["b"], fixed_power=U)

    def build(T):
        s = disk_sample(R, h)
        s.set_boundary_temperature(harmonic_multipole(n, T, R))
        return NonlocalModel(s)

    # 加热增强时卫星环外移 (r_s ∝ Theta), 细网格上结构变化被完全分辨, 所以后半段用较小的步长 + 割线预测
    scan = param_continuation(build, 0.0, M0, s0, [Th * f for f in (1 / 16, 1 / 8, 1 / 4, 3 / 8, 1 / 2, 5 / 8,
                                                                     3 / 4, 7 / 8, 1.0)],
                              lambda model, w0, prev: newton.solve_robust(model, w0, prev["b"], fixed_power=U,
                                                                          maxit=25),
                              secant=True)
    if not scan or scan[-1][0] != Th:
        raise SystemExit(f"未能得到 C_{n} 加热下的定态, 请减小 --Theta")
    _, M3, st3 = scan[-1]
    row = dg.mode_summary(M3, st3, nsym=n, m=m)
    I3 = np.abs(st3["w"]) ** 2
    p_out = float(I3[np.hypot(M3.xo, M3.yo) > 0.5 * R].sum() / I3.sum())
    print(f"start state: b={st3['b']:.4f}, central (small circle) {row['W_small']:+.0f}, "
          f"satellites {row['n_satellites']} at r={row['satellite_r_mean']:.3f}, "
          f"power fraction at r > R/2: {p_out:.1e}, outer singularities {row['n_outer_singularities']}")
    if p_out > 0.01 or row["n_outer_singularities"] > 0:
        print("WARNING: the start state is hybridized with the boundary hot-spot modes (Theta too large for this U);"
              " results are not a clean demonstration - use a smaller --Theta or a larger --U")
    if row["satellite_r_mean"] < 2 * h:
        print(f"WARNING: satellite ring r={row['satellite_r_mean']:.3f} < 2h: not resolved, use a larger --Theta "
              "or a finer grid")
    plotting.plot_mode(M3, st3, os.path.join(out, f"{pre}start_state.png"),
                       title=f"m={m}, C{n} heating Theta={Th:g}, U={U:g}")

    Hc = pattern_field(R, h, n, Th, 0.0)
    Hs = pattern_field(R, h, n, Th, np.pi / (2 * n))     # 旋转 pi/(2n): cos n theta -> sin n theta
    H4 = pattern_field(R, h, 4, Th, 0.0)
    ring = dg.ring_radius(dg.field_grid(M3, st3["w"]), M3.grid, (0, 0), rmax=8.0)

    def rotating(Om):
        # cos(n(theta - Om z)) = cos(n Om z) cos n theta + sin(n Om z) sin n theta, 而 Hs 对应 sin n theta
        return lambda z: np.cos(n * Om * z) * Hc + np.sin(n * Om * z) * Hs

    cases = {"switch_off": lambda z: Hc if z < zs else 0.0 * Hc}
    if n == 3:
        cases["switch_C3_to_C4"] = lambda z: Hc if z < zs else H4
    cases["rotate_slow"] = rotating(0.005 if not args.quick else 0.02)
    cases["rotate_fast"] = rotating(0.05 if not args.quick else 0.2)
    summary = []
    snap = list(np.linspace(0, zmax, 6))
    for name, f in cases.items():
        res = propagate.propagate(M0, st3["w"], zmax, dz=args.dz, diag_every=25, snapshot_zs=snap, nbc_func=f,
                                  w_ref=st3["w"], nsym=n, track_singularities=True, verbose=False)
        plotting.plot_snapshots(M0, res["snapshots"], os.path.join(out, f"{pre}{name}_snapshots.png"),
                                zoom=2 * ring + 1, title=name, show_phase=True)
        tr = satellite_tracks(res, 1.5 * ring)
        io_utils.write_csv(os.path.join(out, f"{pre}{name}_tracks.csv"), tr)
        # C_n 对称性破缺 (非 mu 权重) 与偏离: 区分 "图样驱动的拓扑动力学" 与 "定态不稳定性的放大"
        mu0 = theory.irrep_label(m, n)
        nonmu = [float(sum(w for k, w in irr.items() if k != mu0)) for irr in res["irreps"]]
        io_utils.write_csv(os.path.join(out, f"{pre}{name}_timeseries.csv"),
                           [dict(z=z, nonmu_weight=wn, deviation=d, central=c)
                            for z, wn, d, c in zip(res["z"], nonmu, res["deviation"], res["central"])])
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
        fig.savefig(os.path.join(out, f"{pre}{name}_tracks.png"), dpi=130)
        plt.close(fig)
        qf = dg.field_grid(M0, res["final"])
        nm = np.array(nonmu)
        summary.append(dict(case=name, n=n, m=m, Theta=Th, U=U, zmax=zmax,
                            central_final=round(dg.winding_on_circle(qf, M0.grid, (0, 0), 0.5 * h), 2),
                            power_drift=float(np.ptp(res["power"]) / res["power"][0]),
                            nonmu_initial=float(nm[0]), nonmu_max=float(nm.max()),
                            xi_nonmu_1e3=float(res["z"][np.argmax(nm > 1e-3)]) if (nm > 1e-3).any() else ""))
        print(f"  {name}: final central charge {summary[-1]['central_final']:+.2f}, non-mu weight "
              f"{nm[0]:.1e} -> {nm.max():.1e}  [{res['time']:.0f}s]")
    io_utils.write_csv(os.path.join(out, f"{pre}summary.csv"), summary)
    print("\nE5b done.")


if __name__ == "__main__":
    main()
