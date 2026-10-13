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

2026-10-13 Omega 扫描 (第四轮): 第三轮中双涡旋在 Omega=0.005 时绝热跟随, Omega=0.05 时 Theta=0.5 跟随、Theta=1 在
xi≈164 失锁; 关断后的自由共转角速度 = 圆盘 L=-2 模频率的一半 (0.126)。--omegas 列出若干角速度 (可为负 = 反向旋转),
每个 Omega 一个情形 rotate_Om<+-Omega>; --ramp 让角速度在 [0, ramp] 内线性升到 Omega (避免突然启动激发的振荡,
测的是稳态旋转的临界值); --cases 只跑列出的基本情形 (switch_off/switch_C3_to_C4/rotate_slow/rotate_fast, 给 none 则不跑)。
每个情形另写 <前缀><情形>_row.json (含跟随指标), --table_only 由全部 row.json 重建 control_runs.csv (拆分运行时的汇总表),
并用当前定义由轨迹/时间序列 CSV 重算跟随指标。2026-10-13 (第四轮分析) 改为稳健的定义: 锁定 = 斜坡之后卫星取向转角 / 图样转角
rho ≈ 1 且核心区 (r < 1.5) 奇点集合不变、无相位滑移 (展开的滞后 > pi/n); 原 xi_lock_lost (任何奇点个数变化) 改名 xi_any_change,
它常被外围暗区短暂进出的反涡旋对触发, 不代表失锁 (见 following_metrics)。

运行:
  python e5b_dynamic_control.py --quick    (约 5-10 分钟)
  python e5b_dynamic_control.py            (约 1-2 小时)
  python e5b_dynamic_control.py --hq --n 2 --Theta 0.5 --job e5b_n2_T0.5     (四极加热双涡旋)
  python e5b_dynamic_control.py --hq --n 2 --Theta 1 --cases none --omegas 0.03 --ramp 100 --zmax 600 \
         --tag hq_omega --job n2_T1_Om+0.03                                    (Omega 扫描的一个任务)
  python e5b_dynamic_control.py --hq --tag hq_omega --table_only               (汇总 control_runs.csv)
"""

import glob
import json
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


def pattern_angle(Om, ramp):
    """图样转角 theta0(xi): 角速度在 [0, ramp] 内由 0 线性升到 Omega (ramp <= 0: 突然启动, theta0 = Omega xi)。"""
    if ramp <= 0:
        return lambda z: Om * z
    return lambda z: Om * z * z / (2 * ramp) if z < ramp else Om * (z - ramp / 2)


R_CORE = 1.5          # 核心区半径 (U=15、R=10 时约 0.65 倍亮环半径): 外围暗区的相位噪声伪奇点不算入核心结构
R_RES = 0.15          # 取向只用最内层卫星都在 1.5h 之外的帧 (h=0.1; 更近时位置不可信, §5-2)


def following_metrics(tr, n, angle_fn=None, z_from=None, z_avg=0.0, nonmu=None, r_core=R_CORE, r_res=R_RES):
    """
    由奇点轨迹 (satellite_tracks 的行) 计算图样跟随指标, 对任意 n 适用 (2026-10-13 第四轮分析后改为稳健的定义):
      xi_any_change   r < 1.5 倍亮环半径内的奇点集合 (+1/-1 个数) 第一次改变的 xi —— 原 xi_lock_lost; 第四轮发现它常被
                      外围暗区短暂进出的反涡旋对触发 (双涡旋被拉伸时), 不代表失锁, 只作参考;
      xi_core_change  核心区 r < r_core 内的奇点集合第一次改变的 xi (核心结构改变);
      取向 Phi = arg(sum_k exp(i n phi_k)) / n, 只用最内层 n_plus0 个 +1 卫星, 且只用它们都在 r_res 之外的帧, 展开 (不折回);
      xi_slip         图样旋转时, 展开的滞后 |Phi - Phi0 - theta0(xi)| 第一次超过 pi/n (相位滑移) 的 xi;
      rho             斜坡之后 (xi >= z_avg) 到核心改变为止, 卫星取向转过的角度 / 图样转过的角度 (锁定 = 1);
      locked          rho 与 1 相差 < 2%, 且无核心改变、无相位滑移;
      lag_mean/lag_max 展开的滞后 (到核心改变为止);  r2_ratio  斜坡之后到核心改变为止 r_sat^2 的平均 / 起点的值;
      corotation_rate z_from 之后 (关断情形) 取向的平均角速度 = 自由共转角速度;
      xi_c2           nonmu (C_n 非 mu 权重序列 (z, w)) 第一次超过 1e-3 的 xi (对称性破缺)。
    受迫振子图像 (theory.md §5.1): 双涡旋 = 被加热静态驱动的 L=-2 分裂模; 图样以 Omega 旋转时 r2_ratio = omega_rot/(omega_rot - 2 Omega)。
    """
    by_z = {}
    for t in tr:
        by_z.setdefault(t["z"], []).append(t)
    zs = sorted(by_z)

    def counts(rows, rmax=np.inf):
        return (sum(1 for t in rows if t["charge"] > 0 and t["r"] < rmax),
                sum(1 for t in rows if t["charge"] < 0 and t["r"] < rmax))

    c0 = counts(by_z[zs[0]]) if zs else (0, 0)
    k0 = counts(by_z[zs[0]], r_core) if zs else (0, 0)
    n_sat = max(c0[0], 1)
    any_change = next((z for z in zs if counts(by_z[z]) != c0), "")
    core_change = next((z for z in zs if counts(by_z[z], r_core) != k0), "")
    end = core_change if core_change != "" else np.inf
    z_ok, phi, rsat = [], [], []
    for z in zs:
        if z >= end:
            break
        sat = sorted((t for t in by_z[z] if t["charge"] > 0), key=lambda t: t["r"])[:n_sat]
        if len(sat) < n_sat:
            continue
        rsat.append((z, float(np.mean([t["r"] for t in sat]))))
        if min(t["r"] for t in sat) > r_res:
            z_ok.append(z)
            phi.append(np.angle(sum(np.exp(1j * n * t["angle"]) for t in sat)))
    out = dict(n_plus0=c0[0], n_minus0=c0[1], xi_any_change=any_change, xi_core_change=core_change, xi_slip="",
               rho="", locked="", lag_mean="", lag_max="", corotation_rate="", r2_ratio="",
               r_sat_min=min(r for _, r in rsat) if rsat else "", r_sat_max=max(r for _, r in rsat) if rsat else "",
               xi_c2="")
    if nonmu is not None:
        zc = [z for z, w in nonmu if w > 1e-3]
        out["xi_c2"] = zc[0] if zc else ""
    if rsat:
        zr = np.array([z for z, _ in rsat])
        r2 = np.array([r for _, r in rsat]) ** 2
        if (zr >= z_avg).sum() > 2 and r2[0] > 0:
            out["r2_ratio"] = float(np.mean(r2[zr >= z_avg]) / r2[0])
    if len(z_ok) < 2:
        return out
    z_ok = np.array(z_ok)
    Phi = np.unwrap(np.array(phi)) / n          # phi 已是 n 倍取向 (arg sum exp(i n phi_k)), 展开后除以 n
    if angle_fn is not None:
        th = np.array([angle_fn(z) for z in z_ok])
        lag = Phi - Phi[0] - th
        slip = np.abs(lag) > np.pi / n
        out.update(lag_mean=float(np.mean(lag)), lag_max=float(np.max(np.abs(lag))),
                   xi_slip=float(z_ok[np.argmax(slip)]) if slip.any() else "")
        post = z_ok >= z_avg
        if post.sum() > 10 and abs(th[post][-1] - th[post][0]) > 1e-9:
            out["rho"] = float((Phi[post][-1] - Phi[post][0]) / (th[post][-1] - th[post][0]))
        out["locked"] = bool(out["rho"] != "" and abs(out["rho"] - 1) < 0.02 and core_change == "" and out["xi_slip"] == "")
    if z_from is not None:
        sel = z_ok >= z_from
        if sel.sum() > 2:
            out["corotation_rate"] = float(np.polyfit(z_ok[sel], Phi[sel], 1)[0])
    return out


def assemble(out):
    """汇总全部 control*_row.json; 若对应的轨迹/时间序列 CSV 存在, 用当前的 following_metrics 重算跟随指标 (旧运行也得到新定义)。"""
    rows = []
    for f in sorted(glob.glob(os.path.join(out, "control*_row.json"))):
        with open(f, encoding="utf-8") as fh:
            row = json.load(fh)
        stem = f[:-len("_row.json")]
        ftr, fts = stem + "_tracks.csv", stem + "_timeseries.csv"
        if os.path.exists(ftr):
            tr = [dict(z=float(t["z"]), x=float(t["x"]), y=float(t["y"]), r=float(t["r"]), angle=float(t["angle"]),
                       charge=float(t["charge"])) for t in io_utils.read_csv(ftr)]
            nm = ([(float(t["z"]), float(t["nonmu_weight"])) for t in io_utils.read_csv(fts)]
                  if os.path.exists(fts) else None)
            Om, ramp = float(row.get("omega", 0) or 0), float(row.get("ramp", 0) or 0)
            case = row.get("case", "")
            ang = pattern_angle(Om, ramp) if case.startswith("rotate") else None
            z_from = float(row["zmax"]) / 4 + 50.0 if case == "switch_off" else None
            for k in ("xi_lock_lost",):
                row.pop(k, None)
            row.update(following_metrics(tr, int(row["n"]), angle_fn=ang, z_from=z_from, z_avg=ramp, nonmu=nm))
        rows.append(row)
    rows.sort(key=lambda r: (r.get("n", 0), r.get("Theta", 0), r.get("omega", 0) or 0, r.get("case", "")))
    io_utils.write_csv(os.path.join(out, "control_runs.csv"), rows)
    print(f"control_runs.csv: {len(rows)} runs")


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
    p.add_argument("--cases", nargs="+", default=None,
                   help="只跑这些基本情形 (switch_off switch_C3_to_C4 rotate_slow rotate_fast; none = 一个也不跑); 默认全部")
    p.add_argument("--omegas", nargs="+", type=float, default=None,
                   help="另加以这些角速度旋转的情形 rotate_Om<+-Omega> (负值 = 反向旋转); 用于临界角速度扫描")
    p.add_argument("--ramp", type=float, default=0.0,
                   help="--omegas 情形的角速度在 [0, ramp] 内线性升到 Omega (默认 0: 突然启动, 与 rotate_slow/fast 相同)")
    p.add_argument("--table_only", action="store_true", help="不计算, 只由 control*_row.json 重建 control_runs.csv")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    if args.table_only:
        assemble(out)
        return
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
    row_start = row = dg.mode_summary(M3, st3, nsym=n, m=m)
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

    def rotating(ang):
        # cos(n(theta - a)) = cos(n a) cos n theta + sin(n a) sin n theta, 而 Hs 对应 sin n theta; a = theta0(xi)
        return lambda z: np.cos(n * ang(z)) * Hc + np.sin(n * ang(z)) * Hs

    # 每个情形: (图样函数, 图样转角 theta0(xi) 或 None, 角速度, 启动斜坡长度)
    om_slow, om_fast = (0.005, 0.05) if not args.quick else (0.02, 0.2)
    cases = {"switch_off": (lambda z: Hc if z < zs else 0.0 * Hc, None, 0.0, 0.0)}
    if n == 3:
        cases["switch_C3_to_C4"] = (lambda z: Hc if z < zs else H4, None, 0.0, 0.0)
    for name, Om in (("rotate_slow", om_slow), ("rotate_fast", om_fast)):
        ang = pattern_angle(Om, 0.0)
        cases[name] = (rotating(ang), ang, Om, 0.0)
    if args.cases is not None:
        keep = [c for c in args.cases if c != "none"]
        bad = [c for c in keep if c not in cases]
        if bad:
            raise SystemExit(f"unknown --cases {bad}; available: {list(cases)} or none")
        cases = {k: v for k, v in cases.items() if k in keep}
    for Om in (args.omegas or []):
        ang = pattern_angle(Om, args.ramp)
        cases[f"rotate_Om{Om:+g}"] = (rotating(ang), ang, Om, args.ramp)
    full_set = args.cases is None and not args.omegas
    summary = []
    snap = list(np.linspace(0, zmax, 6))
    for name, (f, ang, Om, ramp) in cases.items():
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
        fm = following_metrics(tr, n, angle_fn=ang, z_from=(zs + 50.0) if name == "switch_off" else None, z_avg=ramp,
                               nonmu=list(zip(res["z"], nonmu)))
        row = dict(summary[-1], omega=Om, ramp=ramp, h=h, R=R, b_start=st3["b"],
                   satellite_r_start=row_start["satellite_r_mean"], **fm)
        with open(os.path.join(out, f"{pre}{name}_row.json"), "w", encoding="utf-8") as fh:
            json.dump(row, fh, indent=1, default=io_utils.json_default)
        print(f"  {name}: final central charge {summary[-1]['central_final']:+.2f}, non-mu weight "
              f"{nm[0]:.1e} -> {nm.max():.1e}; locked={fm['locked']}, rho={fm['rho'] if fm['rho'] == '' else round(fm['rho'], 4)}, "
              f"core change at xi={fm['xi_core_change'] or 'never'}, slip at xi={fm['xi_slip'] or 'never'}, "
              f"co-rotation {fm['corotation_rate'] if fm['corotation_rate'] == '' else round(fm['corotation_rate'], 4)}"
              f"  [{res['time']:.0f}s]", flush=True)
    if full_set:              # 完整运行时保持原来的汇总文件; 拆分运行用 --table_only 汇总 control_runs.csv
        io_utils.write_csv(os.path.join(out, f"{pre}summary.csv"), summary)
    print("\nE5b done.")


if __name__ == "__main__":
    main()
