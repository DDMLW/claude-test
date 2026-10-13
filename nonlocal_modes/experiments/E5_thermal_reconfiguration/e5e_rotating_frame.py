"""
E5e —— 旋转坐标系中的双涡旋: 与旋转加热图样同步转动的定态 (沿 Omega 延拓) 及其线性稳定性
=========================================================================================
目的 (第五轮; theory.md §5.1、Progress.md C10): E5b 的 Omega 扫描 (斜坡启动 + 600 xi 传播) 只给出跟随相图上的离散点,
失锁距离还依赖发射协议与网格。这里直接在旋转坐标系 (model.rotation = Omega, 见 nlmodes/model.py) 中求
“与图样同步转动”的定态 w(Omega) 并算它的 Bogoliubov 谱:
  * 从 Omega=0 (静态加热的双涡旋) 连续延拓得到的定态存在且稳定  <=>  缓慢加速的图样可以被锁定跟随;
    延拓停止 (折返) 或出现增长率 > --growth_tol 的模  ->  理论上的跟随边界 Omega_c(Theta) (与 E5b 的斜坡传播对照);
  * 定态两个 +1 涡旋的半间距给出 r_s^2(Omega)/r_s^2(0), 直接检验受迫分裂模的 omega_rot/(omega_rot - 2 Omega);
  * 旋转坐标系中角向通道 L 的频率平移 L Omega: 分裂模 (L=-2) 在 Omega = omega_{-2}/2、L=-4 模在 omega_{-4}/4 处过零
    (一阶 / 二阶共振), 谱图上可以直接看到哪一支模的碰撞造成失稳 (包括反向旋转时的失稳)。
做法: 半径 R 的圆盘, 边界温度 Theta cos(n theta) (与 E5b/E5d 相同), 固定能流 U; 先沿 Theta 延拓 (Omega=0), 再沿 Omega
正、反两个方向延拓 (割线预测, 失败时步长减半, 连续失败即判为折返); 每隔 --spec_every 个 Omega 点求 Im(lambda) ∈ [0, --wmax]
的低频本征值 (3 个 shift x 16 个, 与 E5d 的 11 x 10 结果相同、快 3 倍), 并在稳定性改变的相邻两点之间补算全部中间点
(边界分辨率 = Omega 步长)。
共振与双稳: 某个被加热驱动的模 (kappa=0) 的旋转系频率 omega_L - L Omega 过零时受迫响应发散; 非线性使响应曲线弯成 S 形
(Duffing 型), 同相支在共振前折返 (hq: Theta=0.5 在 Omega≈0.021, 线性预言的二阶共振 omega_{-4}/4 = 0.028)。延拓失败时在
--max_jump 与其一半处由远到近直接求解 (初值取本段起点与最后一个解), 找到共振另一侧 (反相) 的解后先往回延拓定出这一段的边缘
(若与已有的点重合则说明没有真正的折返), 再往前延拓。各段记在 _branch.csv 的 segment 列; 两段在 Omega 上重叠 = 两支同步定态并存
(两支都稳定处为双稳), 不重叠 = 间隙 (汇总表 branch_overlaps / resonance_gaps; 对应 E5b 中 Theta=0.5、Omega=0.03 失锁而 0.04-0.06 锁定)。

输出 (results/E5_thermal_reconfiguration/<tag>/):
  e5e_n<n>_m<m>_U<U>_T<Theta>_branch.csv  每个 Omega: 所属分段、b、两个 +1 涡旋的半间距 r_s 与 r2_ratio、取向滞后 (相对 Omega=0;
                                          无耗散的受迫响应只有同相 / 反相, 共振以下应为 0)、亮环处谐波比 c1 = c_{m-n}/c_m、
                                          c2 = c_{m-2n}/c_m (一阶 / 二阶受迫响应, 跨过共振时反号)、是否算谱、max growth、
                                          最不稳定模的 Im lambda / L / kappa / Krein
  e5e_n<n>_m<m>_U<U>_T<Theta>_modes.csv   全部低频本征值 (Omega、segment、lam_re、lam_im、L、kappa、krein、w_u)
  e5e_n<n>_m<m>_U<U>_T<Theta>.png         (a) 低频谱 vs Omega  (b) 增长率  (c) r2_ratio 与 omega_rot/(omega_rot - 2 Omega)
  e5e_n<n>_m<m>_U<U>_summary.csv          每个 Theta: 锁定区间 [Omega_lock_min, Omega_lock_max] (从 0 出发, 同相支存在且稳定)、
                                          延拓到达的范围、稳定 / 失稳区间 (lo..hi[段号])、共振间隙与双稳重叠、r2_ratio 对理论的最大偏差
                                          (--table_only 由各 Theta 的 _branch.csv 重建)
运行:
  python e5e_rotating_frame.py --quick                                  (约 3 分钟; h=0.4 只检查流程)
  python e5e_rotating_frame.py --hq --thetas 0.5 --job T0.5              (单核约 2-4 h; 超算按 Theta 拆分, 见 hpc/make_jobs.py)
  python e5e_rotating_frame.py --hq --table_only --job summary          (汇总, 几秒)
"""

import glob
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, heated_disk, initial_solution, param_continuation, pick, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, stability  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from e5d_twin_modes import L_COLORS  # noqa: E402

R_PAIR = 2.0     # 只在 r < 2 (亮环半径约 2.3 以内、光束区内) 找双涡旋的 +1 奇点, 外围暗区的伪奇点不计
T_START = time.time()


def elapsed():
    return f"[{(time.time() - T_START) / 60:.1f} min]"


def pair_geometry(M, w, n):
    """最内层 n 个 +1 奇点: 平均半径 r_s 与 n 重取向 arg(sum e^{i n phi_k})/n。找不到时返回 nan。
    中心荷为 0 (双涡旋) 时不排除中心小圆 (否则 quick 网格上 r_s ≈ 1.25h 的卫星会被漏掉); 否则按惯例排除 1.25h (numerics.md §6)。"""
    g = M.grid
    q = M.to_grid(np.asarray(w, complex))
    reg = dg.beam_region(q, g, rel_threshold=1e-2) & (g.R < R_PAIR)
    central = int(round(dg.winding_on_circle(q, g, (0.0, 0.0), 0.5 * g.h)))
    excl = 1.25 * g.h if central != 0 else 0.0
    sings = [s for s in dg.find_singularities(q, g, reg, exclude_center=(0.0, 0.0), exclude_radius=excl)
             if s[2] > 0]
    sings = sorted(sings, key=lambda s: np.hypot(s[0], s[1]))[:n]
    if len(sings) < n:
        return float("nan"), float("nan")
    r = float(np.mean([np.hypot(x, y) for x, y, _ in sings]))
    axis = float(np.angle(sum(np.exp(1j * n * np.arctan2(y, x)) for x, y, _ in sings)) / n)
    return r, axis


def harmonic_ratios(M, w, r, m, n):
    """亮环半径 r 处的角向谐波比 c_{m-n}/c_m 与 c_{m-2n}/c_m (复数): 分别是被加热一阶驱动的 L=-n 响应与二阶的 L=-2n 响应;
    跨过相应的共振 (旋转系频率过零) 时发散并反号, 用来识别延拓中“一步跨过”的窄共振。"""
    q = M.to_grid(np.asarray(w, complex))
    ls, c = dg.angular_spectrum(q, M.grid, (0.0, 0.0), r, lmax=abs(m) + 2 * n + 2)
    cl = dict(zip(ls.tolist(), c))
    return cl[m - n] / cl[m], cl[m - 2 * n] / cl[m]


def low_spectrum(M, st, shifts, kk, m, n):
    """Im(lambda) >= 0 的低频本征值 (去掉 |lambda| < 1e-5 的相位零模), 按角向通道 L、kappa、Krein 符号标记。"""
    rows, seen = [], []
    for sg in shifts:
        for lam, md in stability.eigenpairs_near(M, st, sg, k=kk):
            if lam.imag < -1e-9 or abs(lam) < 1e-5 or any(abs(lam - x) < 1e-6 for x in seen):
                continue
            seen.append(lam)
            ac = stability.mode_angular_content(M, md, m, w=st["w"])
            rows.append(dict(lam_re=float(lam.real), lam_im=float(lam.imag), L=ac["L"],
                             kappa=((ac["L"] + n // 2) % n) - n // 2, krein=stability.krein_sign(md),
                             w_u=round(ac["w_u"], 3)))
    return rows


def omega_branch(Mt, st, U, om_end, dom, n, m, max_jump, min_frac=1.0 / 16):
    """
    从 Omega=0 沿一个方向延拓到 om_end。目标总是下一个 Omega 网格点 (dom 的整数倍), 失败时步长减半 (割线预测);
    步长小于 dom * min_frac 仍失败 = 折返或共振。此时跳过间隙: 在最后一个解之后 max_jump、max_jump/2 处
    (由远到近, 远处更可能已在共振另一侧) 直接求解, 初值先取本段起点的解、再取最后一个解; 收敛且仍是双涡旋
    (r < 2 内找到 n 个 +1 涡旋) 就开始新的一段: 先往回延拓 (步长 min(dom/2, 0.005), 最多越过近侧折返点 max_jump, 用来定出间隙另一侧的
    边缘; 两段重叠 = 两支同步定态并存; 若与本段已有的点重合, 说明之前的失败只是数值上的, 不另起一段), 再往前延拓。
    返回 [(Omega, model, sol, segment)] (含步长减半时的中间点);
    segment = 0 为含 Omega=0 的一段, 正向依次 +1, +2 ...、反向 -1, -2 ...
    """
    def solve(model, w0, b0, robust=False):
        """延拓的每一步只用 Armijo 牛顿、最多 12 次迭代 (好初值时 <= 6 次收敛; 折返附近快速失败再减半步长);
        跳跃时初值较差, 依次试三种全局化策略 (每种最多 10 次)。"""
        if robust:
            return newton.solve_robust(model, w0, b0, fixed_power=U, maxit=10)
        return newton.solve(model, w0, b0, fixed_power=U, maxit=12)

    def b_guess(frm, v):
        """b 随 Omega 近似平移 m Omega (圆盘中的纯涡旋严格如此, 见单元测试)。"""
        return frm[2]["b"] + m * (v - frm[0])

    def march(start, d, stop, seg, max_step=dom, known=None):
        """从 start 出发沿 d (= ±1) 方向走 Omega 网格点直到 stop (Omega 值) 或失败; 返回 (点列表, 最后一个收敛的状态,
        是否与 known ({Omega: b}, 另一段已有的点) 重合)。重合说明两段是同一支解 (之前的失败只是数值上的, 不是折返)。"""
        pts, cur, prev, step = [], start, None, max_step
        while d * cur[0] < d * stop - 1e-12:
            t = d * cur[0]
            tt = min(t + step, (np.floor(t / dom + 1e-9) + 1) * dom, d * stop)
            vt = round(d * tt, 9)
            Mv = Mt.rotating_copy(vt)
            guess = cur[2]["w"]
            if prev is not None and prev[0] != cur[0]:
                w_prev = newton.align_phase(prev[2]["w"], guess)
                guess = guess + (vt - cur[0]) / (cur[0] - prev[0]) * (guess - w_prev)
            s = solve(Mv, guess, b_guess(cur, vt))
            if s["converged"]:
                prev, cur = cur, (vt, Mv, s)
                pts.append((vt, Mv, s, seg))
                step = min(max_step, 2 * step)
                print(f"   branch Omega={vt:+.5f}: b={s['b']:.6f} ({s['iterations']} Newton its, segment {seg:+d}) "
                      f"{elapsed()}", flush=True)
                kb = (known or {}).get(round(vt, 9))
                if kb is not None and abs(kb - s["b"]) < 1e-7 * max(1.0, abs(kb)):
                    return pts, cur, True
                continue
            step /= 2
            if step < dom * min_frac:
                break
        return pts, cur, False

    sign = 1.0 if om_end > 0 else -1.0
    seg, seg_start = 0, (0.0, Mt, st)
    out, cur, _ = march(seg_start, sign, om_end, seg)
    out = list(out)
    while sign * cur[0] < sign * om_end - 1e-12:
        u0 = sign * cur[0]
        cands = sorted({min((np.ceil((u0 + f * max_jump) / dom - 1e-9)) * dom, abs(om_end))
                        for f in (1.0, 0.5)}, reverse=True)
        nxt = None
        for ug in cands:
            if ug <= u0 + 1e-12:
                continue
            v = round(sign * ug, 9)
            Mj = Mt.rotating_copy(v)
            # 跳跃的初值离解较远: Armijo 线搜索常卡住而纯牛顿收敛 (quick 实测), 所以两种初值都依次试三种全局化策略
            for frm in (seg_start, cur):
                sj = solve(Mj, frm[2]["w"], b_guess(frm, v), robust=True)
                if sj["converged"] and np.isfinite(pair_geometry(Mj, sj["w"], n)[0]):
                    nxt = (v, Mj, sj)
                    break
            if nxt is not None:
                break
        if nxt is None:
            print(f"   co-rotating branch ends at Omega={cur[0]:+g} (fold or resonance wider than {max_jump:g})",
                  flush=True)
            break
        # 往回延拓 (步长 min(dom/2, 0.005), 最多越过近侧折返点 max_jump): 定出另一侧的边缘; 若与本段已有的点重合, 说明根本没有间隙
        known = {round(p[0], 9): p[2]["b"] for p in out + [(0.0, Mt, st, 0)] if p[3] == seg}
        new_seg = seg + int(sign)
        back, _, rejoined = march(nxt, -sign, cur[0] - sign * max_jump, new_seg, max_step=min(dom / 2, 0.005),
                                  known=known)
        if rejoined:
            print(f"   no real gap: the branch beyond Omega={cur[0]:+g} rejoins segment {seg:+d} "
                  f"(Newton had failed there numerically)", flush=True)
            back = [(p[0], p[1], p[2], seg) for p in back[:-1]]
        else:
            seg = new_seg
            print(f"   resonance gap after Omega={cur[0]:+g}: segment {seg:+d} found at {nxt[0]:+g}, extends back to "
                  f"{(back[-1][0] if back else nxt[0]):+g}", flush=True)
        out += back + [(nxt[0], nxt[1], nxt[2], seg)]
        fwd, cur, _ = march(nxt, sign, om_end, seg)
        out += fwd
        seg_start = nxt
    # 去掉重复的 (Omega, 段) 点 (往回延拓可能落在已有的网格点上)
    uniq = {}
    for p in out:
        uniq.setdefault((round(p[0], 9), p[3]), p)
    return list(uniq.values())


def wrap(a, n):
    """把 n 重取向的角度差折回 (-pi/n, pi/n]。"""
    p = np.pi / n
    return float((a + p) % (2 * p) - p)


def theory_r2(Om, w_rot):
    return w_rot / (w_rot - 2 * Om) if w_rot and np.isfinite(w_rot) else float("nan")


def split_mode_frequency(modes0):
    """Omega=0 处转动/分裂模 omega_rot: L=-2、Krein + 的最低稳定支 (与 E5d 的连续追踪起点相同; 强加热时标签可能变成 -4,
    此时退回 kappa=0、Krein + 的最低支)。"""
    st = [r for r in modes0 if r["krein"] == 1 and r["lam_re"] < 1e-5]
    c = sorted(r["lam_im"] for r in st if r["L"] == -2)
    if not c:
        c = sorted(r["lam_im"] for r in st if r["kappa"] == 0)
    return c[0] if c else float("nan")


def run_theta(args, T, Mt, st, name, out):
    n, m, U = args.n, args.m, args.U
    dom = args.dom
    # shift: [0, wmax] 内等距的 ns 个, 每个求最近的 kk 个本征值。云端 h=0.2 检验 (Theta=0.5, Omega=0 与 -0.14):
    # 3 个 shift x 16 与 E5d 式的 11 x 10 找到完全相同的模 (含失稳模 3.08e-2), 快 3 倍 (每个 shift 要一次大 LU)
    ns = pick(args, quick=2, main=3, hq=3)
    shifts = [1j * args.wmax * (2 * j + 1) / (2 * ns) for j in range(ns)]
    kk = pick(args, quick=10, main=16, hq=16)
    every = args.spec_every
    tag = f"{name}_T{T:g}"
    branch = [(0.0, Mt.rotating_copy(0.0), st, 0)]
    for om_end in (args.om_max, args.om_min):
        if om_end != 0:
            branch += omega_branch(Mt, st, U, om_end, dom, n, m, args.max_jump)
    branch.sort(key=lambda t: (t[3], t[0]))
    oms = [t[0] for t in branch]
    segs = sorted({t[3] for t in branch})
    print(f"Theta={T:g}: co-rotating states found for Omega = {min(oms):+g} ... {max(oms):+g} ({len(branch)} points, "
          f"segments " + ", ".join(f"{sg:+d}: {min(t[0] for t in branch if t[3] == sg):+.4f}..."
                                   f"{max(t[0] for t in branch if t[3] == sg):+.4f}" for sg in segs) + ")", flush=True)
    r0, ax0 = pair_geometry(Mt, st["w"], n)
    r_ring = dg.ring_radius(Mt.to_grid(np.asarray(st["w"], complex)), Mt.grid, (0.0, 0.0))
    rows = []
    for Om, M, s, seg in branch:
        r, ax = pair_geometry(M, s["w"], n)
        h1, h2 = harmonic_ratios(M, s["w"], r_ring, m, n)
        rows.append(dict(Theta=T, Omega=Om, segment=seg, b=s["b"], r_s=r,
                         r2_ratio=(r / r0) ** 2 if r0 > 0 else float("nan"),
                         c1_re=h1.real, c1_im=h1.imag, c2_re=h2.real, c2_im=h2.imag,
                         axis_lag=wrap(ax - ax0, n) if np.isfinite(ax) and np.isfinite(ax0) else float("nan"),
                         spectrum=False, max_growth="", unstable_im="", unstable_L="", unstable_kappa="",
                         unstable_krein="", n_modes=""))
    modes = []

    def spectrum_at(i):
        Om, M, s, sg = branch[i]
        ms = low_spectrum(M, s, shifts, kk, m, n)
        for r in ms:
            r.update(Theta=T, Omega=Om, segment=sg)
        modes.extend(ms)
        row = rows[i]
        g = max([r["lam_re"] for r in ms] + [0.0])
        row.update(spectrum=True, max_growth=g, n_modes=len(ms))
        if g > args.growth_tol:
            u = max(ms, key=lambda r: (r["lam_re"], r["lam_im"]))
            row.update(unstable_im=u["lam_im"], unstable_L=u["L"], unstable_kappa=u["kappa"], unstable_krein=u["krein"])
        print(f"   Omega={Om:+.4f}: b={s['b']:.5f}, r_s={row['r_s']:.4f} (r2 ratio {row['r2_ratio']:.3f}), "
              f"lag={row['axis_lag']:+.4f}, {len(ms)} modes, max growth {g:.2e} {elapsed()}", flush=True)

    # 谱只在 Omega 网格点 (dom 的整数倍) 与每一段的两端 (折返点 / 共振间隙两侧) 上算, 逐段处理 (branch 已按 (段, Omega) 排序):
    # 第一遍每隔 every 个网格点; 第二遍在同一段内稳定性改变的相邻两点之间补算全部网格点 (边界分辨率 = dom)
    for sg in segs:
        idx = [i for i in range(len(branch)) if branch[i][3] == sg]
        grid = [i for i in idx if abs(oms[i] / dom - round(oms[i] / dom)) < 1e-6]
        anchor = min(grid, key=lambda i: abs(oms[i])) if grid else idx[0]
        k0 = grid.index(anchor) if grid else 0
        first = sorted({grid[k] for k in range(len(grid)) if (k - k0) % every == 0} | {idx[0], idx[-1]})
        for i in sorted(first, key=lambda i: abs(oms[i])):
            spectrum_at(i)
        done = [i for i in idx if rows[i]["spectrum"]]
        for a, b in zip(done[:-1], done[1:]):
            if (rows[a]["max_growth"] > args.growth_tol) != (rows[b]["max_growth"] > args.growth_tol):
                for i in grid:
                    if a < i < b:
                        spectrum_at(i)
    w_rot = split_mode_frequency([r for r in modes if r["Omega"] == 0.0 and r["segment"] == 0])
    for row in rows:
        row["r2_theory"] = theory_r2(row["Omega"], w_rot)
        row["omega_rot0"] = w_rot
    rows.sort(key=lambda r: (r["Omega"], r["segment"]))
    io_utils.write_csv(os.path.join(out, f"{tag}_branch.csv"), rows)
    io_utils.write_csv(os.path.join(out, f"{tag}_modes.csv"), sorted(modes, key=lambda r: (r["Omega"], r["lam_im"])))
    plot_theta(rows, modes, T, w_rot, os.path.join(out, f"{tag}.png"), args.growth_tol)
    return rows


def plot_theta(rows, modes, T, w_rot, fn, tol):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    ax = axs[0]
    for L in sorted({r["L"] for r in modes}):
        col = L_COLORS.get(L, "0.6")
        for kr, mk in ((1, "o"), (-1, "^"), (0, "s")):
            sel = [r for r in modes if r["L"] == L and r["krein"] == kr]
            if sel:
                ax.plot([r["Omega"] for r in sel], [r["lam_im"] for r in sel], mk, ms=3.5, ls="none", color=col,
                        mfc=col if kr == 1 else "none", mec=col, label=f"L={L:+d}" if kr == 1 else None)
    uns = [r for r in modes if r["lam_re"] > tol]
    ax.plot([r["Omega"] for r in uns], [r["lam_im"] for r in uns], "s", ms=8, mfc="none", mec="#e34948", ls="none",
            label="unstable")
    oo = np.array([r["Omega"] for r in rows])
    if np.isfinite(w_rot):
        ax.plot(oo, np.abs(w_rot - 2 * oo), "--", color="0.35", lw=0.9, label="|omega_rot - 2 Omega|")
    ax.set_ylim(0, None)
    ax.set_xlabel("pattern angular velocity Omega")
    ax.set_ylabel("Im(lambda) in the rotating frame")
    ax.set_title(f"(a) Theta={T:g}: low-lying modes (o Krein +, ^ Krein -)", fontsize=8.5)
    ax.legend(fontsize=6, ncol=3)
    ax = axs[1]
    sp = [r for r in rows if r["spectrum"]]
    g = np.array([max(float(r["max_growth"]), 1e-9) for r in sp])
    ax.semilogy([r["Omega"] for r in sp], g, "o-", ms=3.5, color="#2a78d6")
    ax.axhline(tol, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("pattern angular velocity Omega")
    ax.set_ylabel("max Re(lambda)")
    ax.set_title("(b) growth rate of the co-rotating state (dotted: threshold)", fontsize=8.5)
    ax = axs[2]
    ok = [r for r in rows if np.isfinite(r["r2_ratio"])]
    ax.semilogy([r["Omega"] for r in ok], [r["r2_ratio"] for r in ok], "o", ms=3.5, color="#2a78d6",
                label="co-rotating stationary state")
    if np.isfinite(w_rot):
        o2 = np.linspace(min(oo), min(max(oo), 0.95 * w_rot / 2), 200)
        ax.semilogy(o2, w_rot / (w_rot - 2 * o2), "-", color="0.35", lw=1.0,
                    label=f"omega_rot/(omega_rot - 2 Omega), omega_rot={w_rot:.4f}")
    ax.set_xlabel("pattern angular velocity Omega")
    ax.set_ylabel("r_s^2(Omega) / r_s^2(0)")
    ax.set_title("(c) separation of the two +1 vortices", fontsize=8.5)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(fn, dpi=130)
    plt.close(fig)


def runs_of(sp, tol):
    """sp: 同一段内按 Omega 排好的谱点 [(Omega, growth, row)] -> 状态相同的连续区间 [(状态, [点...])]。"""
    runs = []
    for t in sp:
        st = "unstable" if t[1] > tol else "stable"
        if runs and runs[-1][0] == st:
            runs[-1][1].append(t)
        else:
            runs.append((st, [t]))
    return runs


def summarize(out, name, tol, om_lim=(-0.25, 0.2)):
    """
    各 Theta 的 _branch.csv -> 汇总:
      锁定区间 [Omega_lock_min, Omega_lock_max]: 含 Omega=0 的那一段中、含 0 的稳定区间 (边界取最后一个稳定谱点与第一个
        不稳定谱点的中点; 若一直稳定到段的端点, 取端点, 原因记为 resonance gap / end of branch / scan limit);
      全部稳定区间与失稳区间 (含共振间隙另一侧的分段)、共振间隙、r2_ratio 对 omega_rot/(omega_rot - 2 Omega) 的最大偏差
      (锁定区间内 |Omega| <= 0.1)。
    """
    summ = []
    for f in sorted(glob.glob(os.path.join(out, f"{name}_T*_branch.csv")),
                    key=lambda f: float(re.search(r"_T([0-9.]+)_branch", f).group(1))):
        rows = sorted(io_utils.read_csv(f), key=lambda r: float(r["Omega"]))
        om = np.array([float(r["Omega"]) for r in rows])
        for r in rows:
            r["seg"] = int(float(r.get("segment", 0) or 0))
        segs = sorted({r["seg"] for r in rows})
        rng = {sg: (min(float(r["Omega"]) for r in rows if r["seg"] == sg),
                    max(float(r["Omega"]) for r in rows if r["seg"] == sg)) for sg in segs}
        # 相邻两段之间: 间隙 (没有找到同步定态) 或重叠 (两支同步定态并存; 是否双稳看各自的稳定区间)
        gaps, overlaps = [], []
        for sg in segs:
            nb = sg + (1 if sg >= 0 else -1)
            if sg == 0:
                pairs = [(0, 1), (0, -1)]
            else:
                pairs = [(sg, nb)]
            for a, b in pairs:
                if b not in rng:
                    continue
                if b > 0:
                    lo, hi = rng[a][1], rng[b][0]
                else:
                    lo, hi = rng[b][1], rng[a][0]
                (gaps if lo < hi else overlaps).append((min(lo, hi), max(lo, hi)))
        stable_iv, unstable_iv, lock = [], [], None
        for sg in segs:
            mine = [r for r in rows if r["seg"] == sg]
            o_lo, o_hi = rng[sg]
            sp = [(float(r["Omega"]), float(r["max_growth"]), r) for r in mine if r["spectrum"] == "True"]
            runs = runs_of(sp, tol)
            for k, (st, pts) in enumerate(runs):
                lo = o_lo if k == 0 else 0.5 * (runs[k - 1][1][-1][0] + pts[0][0])
                hi = o_hi if k == len(runs) - 1 else 0.5 * (pts[-1][0] + runs[k + 1][1][0][0])
                if st == "stable":
                    stable_iv.append((lo, hi, sg))
                    if sg == 0 and lo <= 0.0 <= hi:
                        def why(edge, k2, side):
                            if 0 <= k2 < len(runs):
                                u = max(runs[k2][1], key=lambda t: t[1])[2]
                                return f"instability (L={u['unstable_L']}, kappa={u['unstable_kappa']})"
                            lim = om_lim[1] if side > 0 else om_lim[0]
                            if abs(edge - lim) < 1e-6:
                                return "scan limit"
                            return "fold (branch continues beyond the resonance)" if side in rng else "fold (end of branch)"
                        lock = (lo, why(lo, k - 1, -1), hi, why(hi, k + 1, 1))
                else:
                    u = max(pts, key=lambda t: t[1])
                    unstable_iv.append((lo, hi, sg, u[1], u[2]["unstable_L"], u[2]["unstable_kappa"], u[2]["unstable_im"]))
        dev = [float(r["r2_ratio"]) / float(r["r2_theory"]) - 1 for r in rows
               if lock and r["seg"] == 0 and lock[0] <= float(r["Omega"]) <= lock[2] and abs(float(r["Omega"])) <= 0.1
               and r["r2_ratio"] not in ("", "nan") and r["r2_theory"] not in ("", "nan")
               and np.isfinite(float(r["r2_ratio"])) and np.isfinite(float(r["r2_theory"]))]
        summ.append(dict(
            Theta=float(rows[0]["Theta"]), omega_rot0=float(rows[0]["omega_rot0"]),
            Omega_branch_min=float(om.min()), Omega_branch_max=float(om.max()),
            Omega_lock_min=lock[0] if lock else "", lock_min_reason=lock[1] if lock else "",
            Omega_lock_max=lock[2] if lock else "", lock_max_reason=lock[3] if lock else "",
            resonance_gaps="; ".join(f"{a:+.4f}..{b:+.4f}" for a, b in gaps),
            branch_overlaps="; ".join(f"{a:+.4f}..{b:+.4f}" for a, b in overlaps),
            stable_intervals="; ".join(f"{a:+.4f}..{b:+.4f}[{sg:+d}]" for a, b, sg in stable_iv),
            unstable_intervals="; ".join(f"{a:+.4f}..{b:+.4f}[{sg:+d}] (max {g:.2e}, L={L}, kappa={k}, Im={float(im):.3f})"
                                         for a, b, sg, g, L, k, im in unstable_iv),
            n_spectra=sum(1 for r in rows if r["spectrum"] == "True"),
            r2_theory_maxdev_lock=max(np.abs(dev)) if dev else ""))
    io_utils.write_csv(os.path.join(out, f"{name}_summary.csv"), summ)
    return summ


def print_summary(summ):
    for r in summ:
        print(f"Theta={r['Theta']:g}: lock {r['Omega_lock_min']:+.4f} ({r['lock_min_reason']}) ... "
              f"{r['Omega_lock_max']:+.4f} ({r['lock_max_reason']})" if r["Omega_lock_min"] != "" else
              f"Theta={r['Theta']:g}: static state unstable", flush=True)
        print(f"   stable: {r['stable_intervals'] or 'none'}; gaps: {r['resonance_gaps'] or 'none'}; "
              f"two branches: {r['branch_overlaps'] or 'none'}; unstable: {r['unstable_intervals'] or 'none'}", flush=True)


def main():
    p = standard_args(__doc__)
    p.add_argument("--R", type=float, default=10.0)
    p.add_argument("--n", type=int, default=2, help="加热图样阶数 (默认 2: 四极加热双涡旋)")
    p.add_argument("--m", type=int, default=2)
    p.add_argument("--U", type=float, default=15.0, help="固定能流 (与 E5b/E5d 相同)")
    p.add_argument("--thetas", nargs="+", type=float, default=None,
                   help="默认 quick 0.5; main 0.25 0.5 1; hq 0.125-1.1 共 9 个 (超算按 Theta 拆分)")
    p.add_argument("--om_max", type=float, default=0.2)
    p.add_argument("--om_min", type=float, default=-0.25)
    p.add_argument("--dom", type=float, default=None, help="Omega 步长 (默认 quick 0.04 / main 0.01 / hq 0.005)")
    p.add_argument("--spec_every", type=int, default=None,
                   help="每隔几个 Omega 点算一次谱 (默认 quick 1 / main 1 / hq 2); 稳定性改变的两点之间自动补算")
    p.add_argument("--wmax", type=float, default=0.62, help="shift 的最高频率 (Im)")
    p.add_argument("--growth_tol", type=float, default=1e-5, help="判为不稳定的增长率 (hq 的数值地板约 1e-7)")
    p.add_argument("--max_jump", type=float, default=0.03, help="延拓失败时最多跳过的 Omega 间隙 (共振间隙)")
    p.add_argument("--table_only", action="store_true", help="不计算, 由各 Theta 的 _branch.csv 重建 _summary.csv")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    R, h, n, m, U = args.R, args.h, args.n, args.m, args.U
    name = f"e5e_n{n}_m{m}_U{U:g}"
    if args.table_only:
        summ = summarize(out, name, args.growth_tol, (args.om_min, args.om_max))
        print_summary(summ)
        print(f"{name}_summary.csv rebuilt ({len(summ)} Theta)")
        return
    if args.dom is None:
        args.dom = pick(args, quick=0.04, main=0.01, hq=0.005)
    if args.spec_every is None:
        args.spec_every = pick(args, quick=1, main=1, hq=2)
    thetas = sorted(args.thetas if args.thetas is not None else
                    pick(args, quick=[0.5], main=[0.25, 0.5, 1.0],
                         hq=[0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.1]))

    # ---- Omega=0: Theta=0 的定态 (固定能流), 沿 Theta 延拓 (细步长 + 割线预测, 与 E5d 相同)
    M0 = heated_disk(R, h, n, 0.0)
    first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [3.0, 3.5, 2.5], verbose=True)
    s0 = newton.solve_robust(M0, first["w"], first["b"], fixed_power=U)
    targets = sorted({round(x, 6) for T in thetas if T > 0 for x in np.linspace(0.0, T, 9)[1:]})
    scan = param_continuation(lambda v: heated_disk(R, h, n, v), 0.0, M0, s0, targets,
                              lambda model, w0, prev: newton.solve_robust(model, w0, prev["b"], fixed_power=U, maxit=25),
                              secant=True, verbose=False)
    states = {0.0: (M0, s0)}
    for T, Mt, st in scan:
        states[round(T, 6)] = (Mt, st)
    for T in thetas:
        if round(T, 6) not in states:
            print(f"Theta={T:g}: not reached by the Theta continuation", flush=True)
            continue
        Mt, st = states[round(T, 6)]
        run_theta(args, T, Mt, st, name, out)
    print_summary(summarize(out, name, args.growth_tol, (args.om_min, args.om_max)))
    print("\nE5e done.")


if __name__ == "__main__":
    main()
