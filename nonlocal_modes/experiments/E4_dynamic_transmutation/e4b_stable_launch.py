"""
E4b —— 动态嬗变 (二): 稳定区发射、种子扫描、加热样品发射
=========================================================
背景 (Progress.md C9, §5 第 11 条): E4 (hq, xi=400) 中 mu 权重大面积漂移、三角形 m=2 的中心 -1 只在 xi≈85-252
出现、圆盘 m=3 对照破裂。分析结论: 这是 "被种子放大的线性不稳定性" —— E4 的功率 U=10/20 恰好落在 E3 的失稳窗口里
(三角形 m=2: b(U=10)=0.85, b(U=20)=2.80, 窗口上端 b_c≈2.7), 网格 C_4 各向异性/舍入误差作为种子,
非 mu 权重按 exp(2 gamma xi) 增长, 破缺距离 xi_b ≈ ln(1/种子)/(2 gamma)。
本实验用三组发射把 "嬗变" 与 "失稳" 分开, 并给出实验可观测的距离窗口 xi_t < xi < xi_b:
  A) 稳定区发射: 三角形 m=2, U=25、30 (b≈3.9、5.2 > b_c), 另加 1% 噪声 —— 预期: 卫星环形成后长期稳定, mu 守恒;
  B) 种子扫描: 三角形 m=2 U=15 (窗口内, 网格种子 ~1e-5) 与正方形 m=3 U=20 (C_4 与网格相容, 无网格种子),
     输入噪声方差 1e-14 ... 1e-3 —— 预期: 破缺距离与 ln(噪声) 线性, 斜率 -1/(2 gamma), gamma 与 E3 增长率一致
     (三角形中小于网格种子的噪声不改变破缺距离);
  C) 圆盘 m=3 对照 (hq 与 h_scale 1.5 各一次): 破裂距离若与网格无关 -> 舍入误差种子 (L 不是 4 的倍数的模式);
  D) 加热样品发射: C_3 加热圆盘 (R=10, Theta=3) 中 m=2, U=15 (E5b 起点的动态版本; 卫星环 r≈0.43=4h 可直接分辨),
     另加 1% 噪声 —— 其线性稳定性由 E5c 给出。

每个 run:
  1) 求同一 (样品, m, U) 的定态 (先 b=2 处 LG 初值求解, 再沿能流延拓到 U; 加热样品再沿 Theta 延拓),
     得到参考的嬗变核半径 r_s* 与亮环半径; LG 光束宽度取 w0 = r_ring sqrt(2/|m|) (亮环与定态重合, 减小呼吸);
  2) 传播 LG 光束 (可加噪声), 每 25 步记录: C_n 不可约表示权重、中心荷 (多个半径的圆周绕数)、
     中心附近 v 与 m 谐波的振幅比 |c_v/c_m| (r0 = 0.2, 0.5, 1.0)、嬗变核半径 r_s_fit (亚网格谐波拟合)、
     与定态的偏离; 每 100 步检测一次奇点;
  3) 汇总: xi_t50 (r_s_fit 首次达到平台值一半的距离 = 嬗变距离), 平台 r_s 的中位数/离散度 (与 r_s* 比较),
     xi_flip (W(0.5h) 首次等于预言中心荷), xi_b / xi_scan (非 mu 权重首次 > 1e-3 / 1e-2 = 破缺距离),
     xi_10x (首次 > 10 倍初值, 大噪声的 run 用), gamma_fit (非 mu 权重的指数增长率)。

输出 (results/E4_dynamic_transmutation/<tag>/), run 名 e4b_<样品>_m<m>_U<U>[_noise<噪声>]:
  <run>_timeseries.csv        上述全部时间序列
  <run>_timeseries.png        (a) r_s(xi) 与 r_s*; (b) 非 mu 权重 (对数) 与拟合; (c) 中心绕数; (d) 与定态偏离
  <run>_snapshots.png         振幅/相位快照 (标出奇点);  <run>_worldlines.png  奇点到中心的距离随 xi
  <run>_row.json              汇总行 -> e4b_summary.csv
  由 --table_only (或不带 --job 的运行结束时) 重建:
  e4b_summary.csv             每个 run 一行 (并补上 E3/E5c 的线性增长率作对照)
  e4b_seed_scan.csv/png       种子扫描: xi_scan vs 噪声, 由斜率得到 gamma_seed, 与 gamma_fit、E3 增长率比较

运行:
  python e4b_stable_launch.py --quick                                   (约 5 分钟, 三种情形各传播到 xi=20)
  python e4b_stable_launch.py --cases poly3:2:25 --job p3m2U25          (一个 run; 超算按此拆分, hq 约 1-3 h)
  python e4b_stable_launch.py --cases poly4:3:20:1e-6 --job p4m3U20n6
  python e4b_stable_launch.py --table_only
情形写法: 样品:m:U[:噪声方差[:stat]], 样品为 poly<n> (保角半径 --rho)、disk (半径 --rho)、heat<n> (C_n 加热圆盘, 半径 --R,
幅度 --Theta)。末尾加 :stat 表示不发射 LG 光束, 而是发射同功率的定态本身 (加噪声) —— 对照组: 定态线性稳定时非 mu 权重
应保持有界; 若 LG 发射出现慢增长而定态发射没有, 说明慢增长来自 LG 光束的呼吸 (参数共振), 而不是定态的不稳定性
(第二轮 A 组 U=25/30 的 gamma≈0.003, 两套网格一致, 而该处定态线性稳定, 见 Progress.md §6.1)。
"""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (build_sample, finalize_args, heated_disk, initial_solution,  # noqa: E402
                    param_continuation, pick, standard_args)

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, propagate, theory  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# 种子扫描的噪声方差要远小于判据 W_SCAN (噪声的非 mu 权重约为方差的 (n-1)/n), 否则一开始就超过判据;
# 三角形在方形网格上本身有 ~1e-5 量级的种子 (网格的 C_4 各向异性), 比它小的噪声不改变破缺距离 (见 assemble)
ALL_CASES = ["poly3:2:25", "poly3:2:30", "poly3:2:25:1e-2",                                  # A 稳定区
             "poly3:2:15", "poly3:2:15:1e-7", "poly3:2:15:1e-5", "poly3:2:15:1e-3",           # B 三角形种子扫描
             "poly4:3:20:1e-14", "poly4:3:20:1e-11", "poly4:3:20:1e-8", "poly4:3:20:1e-5",    # B 正方形种子扫描
             "disk:3:20",                                                                     # C 圆盘对照
             "heat3:2:15", "heat3:2:15:1e-2"]                                                 # D 加热样品
QUICK_CASES = ["poly3:2:25", "poly4:3:20:1e-8", "heat3:2:15"]
W_BREAK = 1e-3                      # 非 mu 权重超过此值 = 对称性破缺 (与 analysis/a3_e4_growth.py 相同)
W_SCAN = 1e-2                       # 种子扫描用的判据 (仍在线性区: 微扰振幅约为光束的 10%)


def parse_case(c):
    parts = c.split(":")
    kind, m, U = parts[0], int(parts[1]), float(parts[2])
    noise = float(parts[3]) if len(parts) > 3 else 0.0
    launch = parts[4] if len(parts) > 4 else "lg"
    if launch not in ("lg", "stat"):
        raise ValueError(f"launch must be lg or stat: {c}")
    return kind, m, U, noise, launch


def run_name(kind, m, U, noise, launch="lg"):
    return (f"e4b_{kind}_m{m}_U{U:g}" + (f"_noise{noise:g}" if noise > 0 else "")
            + ("_stat" if launch == "stat" else ""))


def build_model(args, kind, h):
    """返回 (模型, 对称阶数 n 或 None, 用于 Theta 延拓的构造函数或 None)。"""
    if kind.startswith("heat"):
        n = int(kind[4:])
        return heated_disk(args.R, h, n, args.Theta), n, (lambda T: heated_disk(args.R, h, n, T))
    sample, nsym = build_sample(kind, h, args.rho)
    return NonlocalModel(sample), nsym, None


def stationary_state(args, kind, M, m, U, theta_builder):
    """同一 (样品, m, U) 的定态: b=2 处 LG 初值 -> 沿能流延拓到 U (加热样品: 先在圆盘中做, 再沿 Theta 延拓)。"""
    M0 = theta_builder(0.0) if theta_builder else M
    first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [2.0, 2.5, 1.5, 3.0], verbose=False)
    if first is None:
        return None
    sols, info = newton.continuation(M0, first["w"], first["power"], U, np.sign(U - first["power"]) * 2.0,
                                     param="P", verbose=False)
    if not sols or abs(sols[-1]["power"] - U) > 1e-6 * U:
        print(f"   stationary state: power continuation failed ({info['stop_reason']})", flush=True)
        return None
    s = sols[-1]
    if theta_builder:
        targets = list(np.linspace(0.0, args.Theta, int(np.ceil(args.Theta / 0.5)) + 1)[1:])
        scan = param_continuation(theta_builder, 0.0, M0, s, targets,
                                  lambda model, w0, prev: newton.solve_robust(model, w0, prev["b"], fixed_power=U,
                                                                              maxit=25),
                                  secant=True, verbose=False)
        if not scan or abs(scan[-1][0] - args.Theta) > 1e-9:
            print("   stationary state: Theta continuation failed", flush=True)
            return None
        s = newton.solve_robust(M, scan[-1][2]["w"], scan[-1][2]["b"], fixed_power=U, maxit=25)
        if not s["converged"]:
            return None
    return s


def make_custom_diag(M, m, v):
    """传播中每次诊断调用: 多半径中心绕数、|c_v/c_m|、主导谐波、r_s_fit。"""
    g, h = M.grid, M.h
    radii_W = {"W_0.5h": 0.5 * h, "W_0.1": 0.1, "W_0.2": 0.2, "W_0.4": 0.4}
    lmax = max(abs(m), abs(v)) + 3

    def diag(z, qg):
        d = {k: dg.winding_on_circle(qg, g, (0.0, 0.0), r) for k, r in radii_W.items()}
        if v != m:
            for r0 in (0.2, 0.5, 1.0):
                ls, c = dg.angular_spectrum(qg, g, (0.0, 0.0), r0, lmax=lmax)
                ls = list(ls)
                d[f"cv_over_cm_r{r0:g}"] = float(abs(c[ls.index(v)]) / max(abs(c[ls.index(m)]), 1e-300))
                d[f"l_dom_r{r0:g}"] = int(ls[int(np.argmax(np.abs(c)))])
            d["r_s_fit"] = dg.estimate_core_radius(qg, g, (0.0, 0.0), m, v, [h, 1.5 * h, 2 * h])
        return d
    return diag


def first_above(z, w, level, w0):
    """W 首次超过 level 的距离; 初始权重 w0 > level/10 时无意义 (一开始就接近判据), 返回 nan。"""
    return float(z[np.argmax(w > level)]) if (w > level).any() and w0 < 0.1 * level else float("nan")


def fit_growth(z, w):
    """
    非 mu 权重 W(xi) 的指数增长率 gamma (W ∝ exp(2 gamma xi)) 与几个破缺距离:
      xi_b    : W 首次 > W_BREAK (1e-3);  xi_scan: W 首次 > W_SCAN (1e-2, 种子扫描用);
      xi_10x  : W 首次 > 10 W0 (加 1% 噪声等大种子的 run 用这个)。
    拟合段: max(1e-14, 3 W0) < W < 1e-2, 且在 W 首次达到 1e-2 之前 (去掉开头噪声弥散的平台和饱和段)。
    """
    w0 = float(w[0])
    zb, zs = first_above(z, w, W_BREAK, w0), first_above(z, w, W_SCAN, w0)
    z10 = float(z[np.argmax(w > 10 * w0)]) if (w > 10 * w0).any() else float("nan")
    sel = (w > max(1e-14, 3 * w0)) & (w < 1e-2)
    if (w >= 1e-2).any():
        sel[int(np.argmax(w >= 1e-2)):] = False
    if sel.sum() < 5:
        return float("nan"), zb, zs, z10
    return float(np.polyfit(z[sel], np.log(w[sel]), 1)[0] / 2), zb, zs, z10


def worldlines(res, ring, fname, title):
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for z, sg in zip(res["z"], res["singularities"]):
        if sg is None:
            continue
        for (x, y, c) in sg:
            r = np.hypot(x, y)
            if r < 2 * ring:
                ax.plot(z, r, ".", ms=2, color="#e34948" if c > 0 else "#2a78d6")
    ax.set_xlabel("xi")
    ax.set_ylabel("distance of unit singularity from centre")
    ax.set_title(f"{title}: red +1, blue -1", fontsize=8)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def timeseries_plot(ts, row, fname, title):
    z = np.array([t["z"] for t in ts])
    fig, axs = plt.subplots(2, 2, figsize=(9.5, 6.2))
    ax = axs[0, 0]
    if "r_s_fit" in ts[0]:
        ax.plot(z, [t["r_s_fit"] for t in ts], "-", color="#2a78d6", lw=1.2, label="r_s(xi), harmonic fit")
        if np.isfinite(row.get("r_s_stationary", np.nan)):
            ax.axhline(row["r_s_stationary"], color="#eb6834", lw=1.2, ls="--", label="stationary r_s*")
        if np.isfinite(row.get("xi_t50", np.nan)):
            ax.axvline(row["xi_t50"], color="0.5", lw=0.8)
        ax.set_ylabel("transmutation core radius r_s")
        ax.legend(fontsize=7)
    ax.set_title("(a) growth of the transmutation core", loc="left", fontsize=8.5)
    ax = axs[0, 1]
    ax.semilogy(z, np.maximum([t["nonmu_weight"] for t in ts], 1e-32), "-", color="#2a78d6", lw=1.2,
                label="non-mu weight")
    if np.isfinite(row.get("gamma_fit", np.nan)) and np.isfinite(row.get("xi_b", np.nan)):
        zz = np.linspace(0, row["xi_b"], 50)
        ax.semilogy(zz, W_BREAK * np.exp(2 * row["gamma_fit"] * (zz - row["xi_b"])), "--", color="#eb6834",
                    lw=1.0, label=f"exp(2 gamma xi), gamma={row['gamma_fit']:.4f}")
    elif np.isfinite(row.get("gamma_fit", np.nan)):
        ax.plot([], [], " ", label=f"fitted gamma={row['gamma_fit']:.4f}")
    ax.axhline(W_BREAK, color="0.5", lw=0.8)
    ax.set_ylabel("weight outside the input irrep mu")
    ax.legend(fontsize=7)
    ax.set_title("(b) C_n symmetry breaking", loc="left", fontsize=8.5)
    ax = axs[1, 0]
    for k, col in (("W_0.5h", "#2a78d6"), ("W_0.1", "#eb6834"), ("W_0.2", "#1baf7a")):
        ax.plot(z, [t[k] for t in ts], "-", color=col, lw=1.0, label=k)
    ax.axhline(row["central_pred"], color="0.5", lw=0.8, ls=":")
    ax.set_ylabel("winding number on a small circle")
    ax.legend(fontsize=7)
    ax.set_title("(c) central charge readout", loc="left", fontsize=8.5)
    ax = axs[1, 1]
    if "deviation" in ts[0]:
        ax.plot(z, [t["deviation"] for t in ts], "-", color="#2a78d6", lw=1.0)
    ax.set_ylabel("|| |q| - |w*| || / ||w*||")
    ax.set_title("(d) distance from the stationary state", loc="left", fontsize=8.5)
    for ax in axs.flat:
        ax.set_xlabel("xi")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def one_run(args, kind, m, U, noise, launch, out):
    name = run_name(kind, m, U, noise, launch)
    h = args.h
    M, nsym, theta_builder = build_model(args, kind, h)
    n_irr = nsym or 4
    v = theory.selection_rule(m, nsym)["central"] if nsym else m
    print(f"\n{name}: {M.describe()}\n   predicted central charge {v}", flush=True)
    s = stationary_state(args, kind, M, m, U, theta_builder)
    row = dict(run=name, sample=kind, n=nsym or 0, m=m, U=U, noise=noise, launch=launch, h=M.h, dz=args.dz,
               Theta=args.Theta if theta_builder else 0.0, central_pred=v)
    if s is not None:
        st = dg.mode_summary(M, s, nsym=nsym, m=m)
        ring = st["ring_radius"]
        row.update(b_stationary=s["b"], r_s_stationary=st.get("r_s_fit", np.nan), ring_stationary=ring,
                   W_small_stationary=st["W_small"])
        print(f"   stationary state: b={s['b']:.4f} ring={ring:.3f} r_s*={row['r_s_stationary']:.4f} "
              f"W(0.5h)={st['W_small']:+.2f}", flush=True)
    else:
        ring = 2.0 * np.sqrt(abs(m) / 2.0)
        row.update(b_stationary=np.nan, r_s_stationary=np.nan, ring_stationary=np.nan)
        print("   WARNING: no stationary reference; LG width 2.0", flush=True)
    w0 = ring * np.sqrt(2.0 / abs(m))
    if launch == "stat":
        if s is None:
            raise SystemExit("stationary launch requested but the stationary state was not found")
        q0 = s["w"].copy()
    else:
        q0 = guesses.lg_vortex(M, m, w0)
    q0 *= np.sqrt(U / M.power(q0))
    if noise > 0:
        q0 = propagate.add_noise(q0, noise, seed=11)
    zmax = args.zmax or pick(args, quick=20.0, main=300.0, hq=600.0)
    snap = list(np.linspace(0.0, zmax, 6))
    res = propagate.propagate(M, q0, zmax, dz=args.dz, diag_every=25, snapshot_zs=snap,
                              w_ref=s["w"] if s is not None else None, nsym=n_irr, track_singularities=True,
                              sing_every=4, custom_diag=make_custom_diag(M, m, v), verbose=False)
    mu0 = theory.irrep_label(m, n_irr)
    ts = []
    for i, z in enumerate(res["z"]):
        irr = res["irreps"][i]
        t = dict(z=z, power=res["power"][i], peak=res["peak"][i], mu_weight=irr[mu0],
                 nonmu_weight=float(sum(w for k, w in irr.items() if k != mu0)), central=res["central"][i])
        if s is not None:
            t["deviation"] = res["deviation"][i]
        sg = res["singularities"][i]
        t["n_singularities"] = "" if sg is None else len(sg)
        t.update(res["custom"][i])
        ts.append(t)
    io_utils.write_csv(os.path.join(out, name + "_timeseries.csv"), ts)
    z = np.array([t["z"] for t in ts])
    wn = np.array([t["nonmu_weight"] for t in ts])
    gam, zb, zs, z10 = fit_growth(z, wn)
    row.update(w0=w0, zmax=zmax, gamma_fit=gam, xi_b=zb, xi_scan=zs, xi_10x=z10, nonmu_initial=float(wn[0]),
               nonmu_max=float(wn.max()),
               mu_weight_drift=float(max(abs(t["mu_weight"] - ts[0]["mu_weight"]) for t in ts)),
               power_drift=float(np.ptp(res["power"]) / res["power"][0]), time=round(res["time"], 1))
    # 嬗变距离与平台: 只用破缺之前的部分 (去掉开头 10%)
    if v != m:
        rs = np.array([t["r_s_fit"] for t in ts], float)
        z_end = zb if np.isfinite(zb) else zmax
        plat = (z >= 0.1 * z_end) & (z <= z_end) & np.isfinite(rs)
        r_plat = float(np.median(rs[plat])) if plat.any() else float("nan")
        row.update(r_s_plateau=r_plat, r_s_plateau_std=float(np.std(rs[plat])) if plat.any() else float("nan"))
        reach = np.isfinite(rs) & (rs >= 0.5 * r_plat)
        row["xi_t50"] = float(z[np.argmax(reach)]) if np.isfinite(r_plat) and reach.any() else float("nan")
        flip = np.array([round(t["W_0.5h"]) == v for t in ts])
        row["xi_flip"] = float(z[np.argmax(flip)]) if flip.any() else float("nan")
    qf = dg.field_grid(M, res["final"])
    row["W_small_final"] = round(dg.winding_on_circle(qf, M.grid, (0, 0), 0.5 * M.h), 2)
    print(f"   [{res['time']:.0f}s] xi_t50={row.get('xi_t50', float('nan')):.1f} "
          f"r_s plateau={row.get('r_s_plateau', float('nan')):.4f} (stationary {row['r_s_stationary']:.4f}) "
          f"xi_b={zb:.1f} xi_10x={z10:.1f} gamma_fit={gam:.4f} nonmu {wn[0]:.1e}->{wn.max():.1e} "
          f"power drift {row['power_drift']:.1e}", flush=True)
    title = (f"{name}: " + (f"stationary state m={m}" if launch == "stat" else f"LG m={m} (w0={w0:.2f})")
             + f" into {kind}, U={U:g}" + (f", noise {noise:g}" if noise else ""))
    timeseries_plot(ts, row, os.path.join(out, name + "_timeseries.png"), title)
    plotting.plot_snapshots(M, res["snapshots"], os.path.join(out, name + "_snapshots.png"), zoom=2 * ring + 1.5,
                            title=title, show_phase=True)
    worldlines(res, ring, os.path.join(out, name + "_worldlines.png"), name)
    with open(os.path.join(out, name + "_row.json"), "w", encoding="utf-8") as fh:
        json.dump(row, fh, indent=1, default=io_utils.json_default)


# ---------------------------------------------------------------- 汇总
def linear_growth(row, tag_root):
    """对照: 同一定态的线性增长率 (多边形/圆盘查 E3 谱 (rho=8), 加热样品查 E5c 的 U_mark 结果)。"""
    b = row.get("b_stationary")
    if b is None or not np.isfinite(b):
        return float("nan")
    if row["sample"].startswith("heat"):
        f = os.path.join(tag_root, "E5_thermal_reconfiguration", "hq",
                         f"e5c_n{row['n']}_m{row['m']}_T{row['Theta']:g}_row.json")
        if os.path.exists(f):
            with open(f, encoding="utf-8") as fh:
                r = json.load(fh)
            if abs(r.get("U_mark", -1) - row["U"]) < 1e-9 and "growth_at_U_mark" in r:
                return float(r["growth_at_U_mark"])
        return float("nan")
    f = os.path.join(tag_root, "E3_stability", "hq", f"{row['sample']}_m{row['m']}_spectrum.csv")
    if not os.path.exists(f):
        return float("nan")
    sr = io_utils.read_csv(f)
    bs = np.array([float(r["b"]) for r in sr])
    gs = np.array([float(r["max_growth"]) for r in sr])
    o = np.argsort(bs)
    return float(np.interp(b, bs[o], gs[o])) if bs.min() <= b <= bs.max() else float("nan")


def assemble(out):
    rows = []
    for f in glob.glob(os.path.join(out, "e4b_*_row.json")):
        with open(f, encoding="utf-8") as fh:
            rows.append(json.load(fh))
    if not rows:
        print("no e4b_*_row.json found")
        return
    root = os.path.join(io_utils.PROJECT_ROOT, "results")
    for r in rows:
        r["linear_growth_E3_or_E5c"] = linear_growth(r, root)
    rows.sort(key=lambda r: (r["sample"], r["m"], r["U"], r["noise"]))
    io_utils.write_csv(os.path.join(out, "e4b_summary.csv"), rows)
    print(f"e4b_summary.csv: {len(rows)} runs")
    # 种子扫描: 同一 (样品, m, U, h) 的 xi_scan (非 mu 权重首次 > 1e-2) vs ln(噪声)
    groups = {}
    for r in rows:
        groups.setdefault((r["sample"], r["m"], r["U"], r["h"], r.get("launch", "lg")), []).append(r)
    scan_rows = []
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    styles = [("#2a78d6", "o"), ("#eb6834", "s"), ("#1baf7a", "^"), ("#eda100", "D")]
    k = 0
    for (kind, m, U, h, launch), rs in sorted(groups.items()):
        base = [r for r in rs if r["noise"] == 0]
        xb0 = base[0].get("xi_scan", float("nan")) if base else float("nan")
        pts = sorted([(r["noise"], r["xi_scan"]) for r in rs
                      if r["noise"] > 0 and np.isfinite(r.get("xi_scan", np.nan))])
        # 比网格/舍入种子还小的噪声不改变破缺距离 (被 "种子下限" 截断), 不参加斜率拟合
        pts_fit = [p for p in pts if not (np.isfinite(xb0) and p[1] > 0.95 * xb0)]
        if len(pts_fit) < 2:
            continue
        x = np.log([p[0] for p in pts_fit])
        slope = np.polyfit(x, [p[1] for p in pts_fit], 1)[0]
        gam_seed = -1.0 / (2 * slope) if slope < 0 else float("nan")
        gfit = np.nanmedian([r["gamma_fit"] for r in rs]) if any(np.isfinite(r["gamma_fit"]) for r in rs) else np.nan
        scan_rows.append(dict(sample=kind, m=m, U=U, h=h, n_seeds=len(pts), n_fit=len(pts_fit),
                              slope_dxib_dlnnoise=slope, gamma_seed=gam_seed, gamma_fit_median=gfit,
                              linear_growth=rs[0]["linear_growth_E3_or_E5c"], xi_b_no_noise=xb0))
        col, mk = styles[k % len(styles)]
        k += 1
        ax.semilogx([p[0] for p in pts], [p[1] for p in pts], mk, ls="-", color=col,
                    label=f"{kind} m={m} U={U:g}: gamma_seed={gam_seed:.3f}")
        if np.isfinite(xb0):
            ax.axhline(xb0, color=col, lw=0.8, ls=":")
    io_utils.write_csv(os.path.join(out, "e4b_seed_scan.csv"), scan_rows)
    if scan_rows:
        ax.set_xlabel("input noise variance (seed)")
        ax.set_ylabel("distance where the non-mu weight exceeds 1e-2")
        ax.set_title("seeded instability: xi = ln(1/seed)/(2 gamma) + const\n(dotted: no added noise)",
                     fontsize=8.5, loc="left")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "e4b_seed_scan.png"), dpi=140)
    plt.close(fig)
    for r in scan_rows:
        print(f"seed scan {r['sample']} m={r['m']} U={r['U']:g}: gamma_seed={r['gamma_seed']:.4f} "
              f"gamma_fit={r['gamma_fit_median']:.4f} linear={r['linear_growth']:.4f}")


def main():
    p = standard_args(__doc__)
    p.add_argument("--cases", nargs="+", default=None, help="样品:m:U[:噪声方差] 列表 (默认全部, quick 三种)")
    p.add_argument("--rho", type=float, default=8.0, help="多边形/圆盘的保角半径 (与 E2-E4 相同)")
    p.add_argument("--R", type=float, default=10.0, help="加热圆盘半径 (与 E5 相同)")
    p.add_argument("--Theta", type=float, default=3.0, help="加热幅度 (heat<n> 情形; 与 E5b 起点相同)")
    p.add_argument("--zmax", type=float, default=None, help="传播距离 (默认 quick 20, main 300, hq 600)")
    p.add_argument("--table_only", action="store_true", help="不计算, 只由 e4b_*_row.json 重建汇总表与种子扫描图")
    args = finalize_args(p, "E4_dynamic_transmutation")
    out = args.res_dir
    if args.table_only:
        assemble(out)
        return
    for c in args.cases or (QUICK_CASES if args.quick else ALL_CASES):
        one_run(args, *parse_case(c), out)
    if not args.job:
        assemble(out)
    print("\nE4b done.")


if __name__ == "__main__":
    main()
