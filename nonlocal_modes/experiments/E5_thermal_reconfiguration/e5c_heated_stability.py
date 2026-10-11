"""
E5c —— 加热态的线性稳定性: (Theta, b) 平面的失稳窗口图
=====================================================
背景 (Progress.md C8, 开放问题 Q6、Q9):
  * E3 (hq) 发现 m=2 在所有样品中都有低功率振荡失稳窗口 (圆盘 b_c≈0.53, rho=8), 三角形把它拉宽到 b_c≈2.7,
    三角形中 m=1 也出现窗口。E3b 的机理假说: C_n 边界把角向通道 L 与 L±n 耦合, 两支 Krein 符号相反的频率
    相遇处 (Hamilton-Hopf 碰撞) 打开新的窗口。
  * 多边形只有离散的 n 和固定的各向异性强度, 无法检验 "耦合强度 -> 窗口" 的定量关系。热边界 Theta cos(n theta)
    提供连续可调的 C_n 耦合: 小 Theta 时耦合 ∝ Theta, Hamilton-Hopf 理论预言新窗口出现在圆盘通道交叉点 b*
    附近, 窗口宽度与最大增长率都 ∝ Theta (双对数斜率 1)。
  * E5b 的起点 (n=3, m=2, Theta=3, U=15, b≈2.27) 是否线性稳定 (Q9)? 若不稳定, E5b 长距离演化要重新解读。

本脚本对每个情形 n:m (C_n 加热, 拓扑荷 m) 与每个 Theta:
  1) 在 b_start 处求圆盘 (Theta=0) 的 m 荷涡旋, 固定 b 沿 Theta 延拓到目标值 (割线预测);
  2) 沿 b 向上、向下延拓得到解族 (b 网格是 db 的整数倍, 不同 Theta、不同网格可逐点比较);
  3) 在 spec_step 的整数倍处计算稳定性谱 (与 E3b 相同的低功率 shift 扫描), 稳定/不稳定交界处二分加密;
     每个点同时记录外围功率比例 P_outer (与边缘热点模式杂化的判据, 见 Progress.md §5 第 9 条);
  4) 若给出 --U_mark (m=2 默认 15, 即 E5b 的能流), 另在固定能流处求解并算谱 —— 直接回答 E5b 起点是否稳定。
  Theta=0 就是半径 R 的圆盘; 由标度不变性 (x -> x/R, b -> b R^2), 它的窗口应等于 E3 圆盘 (rho=8) 的窗口乘 (8/R)^2,
  这是本实验自带的一致性检验 (R=10: 圆盘 m=2 的 b_c ≈ 0.53 x 0.64 ≈ 0.34)。
  Theta 反号等价于把图样转 pi/n, 稳定性不变, 所以只算 Theta >= 0。

输出 (results/E5_thermal_reconfiguration/<tag>/):
  e5c_n<n>_m<m>_T<Theta>_spectrum.csv    每个 b: 增长率、最不稳定本征值、L、kappa、P_outer (含二分加密的点)
  e5c_n<n>_m<m>_T<Theta>_row.json        汇总行 (阈值、窗口、U_mark 处的结果)
  e5c_n<n>_m<m>_T<Theta>_growth.png      增长率 vs b (标出 L)
  e5c_n<n>_m<m>_T<Theta>_mode_b<val>.png 每个窗口最不稳定点的模式图 (及 U_mark 处不稳定时的模式图)
  由 --table_only (或不带 --job 的运行结束时) 从上面的分文件重建:
  e5c_summary.csv         每个 (n, m, Theta) 一行
  e5c_windows.csv         每个失稳窗口一行 (边界、宽度、最大增长率、L/kappa、是否为加热诱导的新窗口、是否杂化)
  e5c_stability_map.png   (Theta, b) 平面稳定性图 + 诱导窗口的最大增长率/宽度 vs Theta (双对数, 斜率 1 参考线)

运行:
  python e5c_heated_stability.py --quick                                    (约 10 分钟: 3:2, Theta = 0, 2, 粗网格)
  python e5c_heated_stability.py --cases 3:2 --thetas 1.5 --job n3_m2_T1.5  (一个 (情形, Theta); 超算按此拆分)
  python e5c_heated_stability.py --table_only                               (拆分运行结束后重建汇总表与总图)
"""

import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (finalize_args, heated_disk, initial_solution, outer_power_fraction,  # noqa: E402
                    param_continuation, pick, standard_args)
from stability_tools import (TOL, bisect_thresholds, low_power_shifts, plot_mode, spectrum_point,  # noqa: E402
                             windows_and_thresholds)

from nlmodes import guesses, io_utils, newton  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

ALL_CASES = ["3:2", "3:1", "4:2"]
B_RANGE = {1: (0.05, 1.5), 2: (0.05, 3.0)}      # R=10 时的 b 范围 (约为 E3 rho=8 范围乘 0.64)
B_START = {1: 1.5, 2: 2.5}                      # 沿 Theta 延拓所在的 b (取范围上端: 光束窄、离边界远, 与边缘热点模式杂化所需的 Theta 最大)
THETA_STEP = 0.5                                # 沿 Theta 延拓的最大步长
# 图中情形的固定配色/标记 (与 analysis/style.py 的调色板同一顺序, 已做色盲校验)
CASE_STYLE = {"3:2": ("#2a78d6", "o"), "3:1": ("#eb6834", "s"), "4:2": ("#1baf7a", "^")}


def theta_list(args):
    return pick(args, quick=[0.0, 2.0], main=[0.0, 0.5, 1.0, 2.0, 3.0],
                hq=[0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])


def run_name(n, m, T):
    return f"e5c_n{n}_m{m}_T{T:g}"


def family_at_theta(args, n, m, T, b_lo, b_hi, db, R):
    """Theta 下的解族: b_start 处求圆盘解 -> 固定 b 沿 Theta 延拓 -> 沿 b 上下延拓。返回 (模型, 解列表, 信息)。"""
    b0 = round(float(np.clip(B_START[m], b_lo, b_hi)) / db) * db
    M0 = heated_disk(R, args.h, n, 0.0)
    s0 = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [b0], verbose=False)
    if s0 is None:
        return None, [], dict(stop="no initial solution at Theta=0")
    M, sT = M0, s0
    if T != 0:
        targets = list(np.linspace(0.0, T, int(np.ceil(T / THETA_STEP)) + 1)[1:])
        scan = param_continuation(lambda v: heated_disk(R, args.h, n, v), 0.0, M0, s0, targets,
                                  lambda model, w0, prev: newton.solve_robust(model, w0, prev["b"], maxit=25),
                                  secant=True, verbose=False)
        if not scan or abs(scan[-1][0] - T) > 1e-9:
            reached = scan[-1][0] if scan else 0.0
            return None, [], dict(stop=f"Theta continuation stopped at {reached:g} (b={b0})")
        M, sT = scan[-1][1], scan[-1][2]
    up, info_up = newton.continuation(M, sT["w"], b0, b_hi, db, verbose=False)
    down, info_down = newton.continuation(M, sT["w"], b0, b_lo, -db, verbose=False)
    sols = sorted(down[1:] + up, key=lambda s: s["b"])
    return M, sols, dict(stop_up=info_up["stop_reason"], stop_down=info_down["stop_reason"])


def on_grid(b, step):
    return abs(b / step - round(b / step)) < 1e-6


def spectrum_at_power(M, sols, U, m, n, shifts, kk, R):
    """固定能流 U 处的解与谱 (由能流最接近的解族成员出发)。"""
    near = min(sols, key=lambda s: abs(s["power"] - U))
    if not (min(s["power"] for s in sols) <= U <= max(s["power"] for s in sols)):
        return None, None, None
    sU = newton.solve_robust(M, near["w"], near["b"], fixed_power=U, maxit=25)
    if not sU["converged"]:
        return None, None, None
    row, mode = spectrum_point(M, sU, m, n, shifts, kk)
    row["P_outer"] = outer_power_fraction(M, sU["w"], R)
    return sU, row, mode


def growth_plot(rows, fname, title):
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.plot([r["b"] for r in rows], [max(r["growth"], 0.0) for r in rows], "o-", ms=3, color="#2a78d6")
    for r in rows:
        if r["unstable"] and r.get("L", "") != "":
            ax.annotate(f"{int(r['L']):+d}", (r["b"], r["growth"]), fontsize=6, xytext=(0, 4),
                        textcoords="offset points", ha="center")
    ax.axhline(TOL, color="0.5", lw=0.8)
    ax.set_xlabel("b")
    ax.set_ylabel("max Re(lambda)  (labels: angular index L)")
    ax.set_title(title, fontsize=7)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- 汇总 (--table_only)
def load_spectrum(f):
    rows = []
    for r in io_utils.read_csv(f):
        rows.append(dict(b=float(r["b"]), U=float(r["U"]), growth=float(r["growth"]), lam_im=float(r["lam_im"]),
                         unstable=r["unstable"] == "True", L=r.get("L", ""), kappa=r.get("kappa", ""),
                         P_outer=float(r["P_outer"]) if r.get("P_outer", "") != "" else np.nan))
    return sorted(rows, key=lambda r: r["b"])


def window_edges(rows, win):
    """窗口边界取最外侧不稳定点与相邻稳定点的中点 (二分加密后误差约 db/16); 碰到扫描范围端点则取端点。"""
    i0, i1 = rows.index(win[0]), rows.index(win[-1])
    lo = 0.5 * (rows[i0 - 1]["b"] + win[0]["b"]) if i0 > 0 else win[0]["b"]
    hi = 0.5 * (rows[i1 + 1]["b"] + win[-1]["b"]) if i1 < len(rows) - 1 else win[-1]["b"]
    return lo, hi, i0 == 0, i1 == len(rows) - 1


def assemble(out):
    summ = []
    for f in glob.glob(os.path.join(out, "e5c_n*_row.json")):
        with open(f, encoding="utf-8") as fh:
            summ.append(json.load(fh))
    if not summ:
        print("no e5c_*_row.json found")
        return
    summ.sort(key=lambda r: (r["n"], r["m"], r["Theta"]))
    io_utils.write_csv(os.path.join(out, "e5c_summary.csv"), summ)
    spectra = {}
    for f in glob.glob(os.path.join(out, "e5c_n*_spectrum.csv")):
        n, m, T = re.findall(r"e5c_n(\d+)_m(\d+)_T([\d.]+)_spectrum", os.path.basename(f))[0]
        spectra[(int(n), int(m), float(T))] = load_spectrum(f)
    wrows = []
    for (n, m, T), rows in sorted(spectra.items()):
        ref = spectra.get((n, m, 0.0))
        ref_wins = []
        if ref is not None:
            ref_wins = [window_edges(ref, w)[:2] for w in windows_and_thresholds(ref)[0]]
        for w in windows_and_thresholds(rows)[0]:
            lo, hi, at_lo, at_hi = window_edges(rows, w)
            top = max(w, key=lambda r: r["growth"])
            induced = "" if ref is None else (T > 0 and not any(lo <= b2 and hi >= b1 for b1, b2 in ref_wins))
            hyb = ""
            if ref is not None:
                pref = {round(r["b"], 6): r["P_outer"] for r in ref}
                hyb = any(r["P_outer"] > 2 * pref[round(r["b"], 6)] + 2e-3 for r in w if round(r["b"], 6) in pref)
            wrows.append(dict(n=n, m=m, Theta=T, b_lo=lo, b_hi=hi, width=hi - lo, touches_range_end=at_lo or at_hi,
                              max_growth=top["growth"], b_at_max=top["b"], lam_im_at_max=top["lam_im"],
                              L_at_max=top["L"], kappa_at_max=top["kappa"], induced=induced, hybrid=hyb))
    io_utils.write_csv(os.path.join(out, "e5c_windows.csv"), wrows)
    print(f"e5c_summary.csv: {len(summ)} runs; e5c_windows.csv: {len(wrows)} windows")
    plot_map(out, summ, spectra, wrows)


def plot_map(out, summ, spectra, wrows):
    cases = [c for c in ALL_CASES if any((n, m) == tuple(map(int, c.split(":"))) for n, m, _ in spectra)]
    if not cases:
        return
    g_all = [r["growth"] for rows in spectra.values() for r in rows if r["unstable"]]
    norm = LogNorm(vmin=max(min(g_all), TOL), vmax=max(g_all)) if g_all else LogNorm(TOL, 1e-1)
    cmap = plt.get_cmap("Blues")
    ncol = max(len(cases), 2)
    fig = plt.figure(figsize=(3.9 * ncol, 7.2), layout="constrained")
    gs = fig.add_gridspec(2, ncol)
    sc, map_axes = None, []
    for i, c in enumerate(cases):
        n, m = map(int, c.split(":"))
        ax = fig.add_subplot(gs[0, i])
        map_axes.append(ax)
        Ts = sorted(T for (nn, mm, T) in spectra if (nn, mm) == (n, m))
        for T in Ts:
            rows = spectra[(n, m, T)]
            st = [r for r in rows if not r["unstable"]]
            un = [r for r in rows if r["unstable"]]
            ax.plot([r["b"] for r in st], [T] * len(st), "o", ms=2.5, color="#c9c8c3", mec="none")
            sc = ax.scatter([r["b"] for r in un], [T] * len(un), c=[r["growth"] for r in un], cmap=cmap, norm=norm,
                            s=22, marker="s", edgecolors="none")
        mark = sorted((r["Theta"], r["b_at_U_mark"]) for r in summ if (r["n"], r["m"]) == (n, m)
                      and r.get("b_at_U_mark", "") not in ("", None))
        if mark:
            ax.plot([b for _, b in mark], [T for T, _ in mark], "--", color="#52514e", lw=1.0)
            T3 = [x for x in mark if abs(x[0] - 3.0) < 1e-9]
            if (n, m) == (3, 2) and T3:
                ax.plot(T3[0][1], 3.0, "*", ms=12, color="#eb6834", mec="#0b0b0b", mew=0.6, zorder=5)
                ax.annotate("E5b start (U=15)", (T3[0][1], 3.0), xytext=(6, -14), textcoords="offset points",
                            fontsize=7, color="#0b0b0b")
            ax.annotate(f"U={summ[0].get('U_mark', 15):g}", (mark[0][1], mark[0][0]), xytext=(4, 4),
                        textcoords="offset points", fontsize=7, color="#52514e")
        ax.set_xlabel("propagation constant b")
        ax.set_ylabel("boundary heating Theta")
        ax.set_title(f"({'abc'[i]}) C{n} heating, m={m}", loc="left", fontsize=9)
        ax.set_yticks(Ts)
        ax.set_ylim(min(Ts) - 0.2, max(Ts) + 0.3)
    if sc is not None:
        cb = fig.colorbar(sc, ax=map_axes, shrink=0.9, pad=0.01)
        cb.set_label("growth rate max Re(lambda)\n(squares: unstable, grey: stable)")
    ind = [w for w in wrows if w["induced"] is True and w["Theta"] > 0]
    for j, (key, lab) in enumerate((("max_growth", "max growth rate in the window"), ("width", "window width in b"))):
        ax = fig.add_subplot(gs[1, j])
        for c in cases:
            n, m = map(int, c.split(":"))
            col, mk = CASE_STYLE.get(c, ("#6250d6", "D"))
            sel = sorted((w for w in ind if (w["n"], w["m"]) == (n, m)), key=lambda w: w["Theta"])
            if sel:
                ax.loglog([w["Theta"] for w in sel], [w[key] for w in sel], mk, ls="none", color=col, ms=6,
                          label=f"C{n} heating, m={m}")
        if ind:
            T_ = np.array([0.2, 3.5])
            y0 = np.median([w[key] / w["Theta"] for w in ind])
            ax.loglog(T_, y0 * T_, "-", color="#52514e", lw=1.0, label="slope 1 (Hamilton-Hopf)")
            ax.legend(loc="upper left", fontsize=7)
        ax.set_xlim(0.2, 3.5)
        ax.set_xticks([0.25, 0.5, 1, 2, 3])
        ax.set_xticklabels(["0.25", "0.5", "1", "2", "3"])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("boundary heating Theta")
        ax.set_ylabel(lab)
        ax.set_title(f"({'de'[j]}) heating-induced windows", loc="left", fontsize=9)
    fig.savefig(os.path.join(out, "e5c_stability_map.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- 主程序
def main():
    p = standard_args(__doc__)
    p.add_argument("--cases", nargs="+", default=None, help="n:m 列表 (C_n 加热, 拓扑荷 m), 默认 3:2 3:1 4:2")
    p.add_argument("--thetas", nargs="+", type=float, default=None, help="边界加热幅度列表 (默认按精度档位)")
    p.add_argument("--R", type=float, default=10.0, help="圆盘半径 (与 E5a/E5b 相同)")
    p.add_argument("--b_range", nargs=2, type=float, default=None, help="只算这个 b 区间 (阈值附近加密复核用)")
    p.add_argument("--U_mark", type=float, default=15.0, help="m=2 时另在此固定能流处算谱 (E5b 起点; <=0 关闭)")
    p.add_argument("--table_only", action="store_true", help="不计算, 只由分文件重建 e5c_summary/windows 与总图")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    if args.table_only:
        assemble(out)
        return
    R = args.R
    cases = [tuple(map(int, c.split(":"))) for c in (args.cases or (["3:2"] if args.quick else ALL_CASES))]
    thetas = args.thetas if args.thetas is not None else theta_list(args)
    db = pick(args, quick=0.2, main=0.1, hq=0.05)
    spec_step = pick(args, quick=0.2, main=0.1, hq=0.1)
    shifts = low_power_shifts(args)
    kk = pick(args, quick=10, main=14, hq=18)
    n_bisect = pick(args, quick=1, main=3, hq=4)
    for n, m in cases:
        b_lo, b_hi = args.b_range or B_RANGE[m]
        if args.quick:
            b_lo = max(b_lo, 0.2)
        for T in thetas:
            name = run_name(n, m, T)
            print(f"\n{name}: disk R={R:g}, T_b = {T:g} cos({n} theta), m={m}, b in [{b_lo}, {b_hi}], db={db}",
                  flush=True)
            M, sols, info = family_at_theta(args, n, m, T, b_lo, b_hi, db, R)
            if not sols:
                print(f"   no family: {info}", flush=True)       # 也写一行汇总 (例如 Theta 延拓因杂化中断)
                with open(os.path.join(out, f"{name}_row.json"), "w", encoding="utf-8") as fh:
                    json.dump(dict(n=n, m=m, Theta=T, R=R, h=args.h, n_points=0, stop=info.get("stop", "")), fh,
                              indent=1, default=io_utils.json_default)
                continue
            print(f"   {M.describe()}\n   family: {len(sols)} solutions, b in [{sols[0]['b']:.3f}, {sols[-1]['b']:.3f}]"
                  f"; up: {info['stop_up']}; down: {info['stop_down']}", flush=True)
            pts = [s for s in sols if on_grid(s["b"], spec_step)]
            rows, modes = [], {}
            for s in pts:
                row, mode = spectrum_point(M, s, m, n, shifts, kk)
                row["P_outer"] = outer_power_fraction(M, s["w"], R)
                rows.append(row)
                modes[round(s["b"], 6)] = (s, mode)
                print(f"   b={s['b']:.3f} U={s['power']:.3f} growth={row['growth']:.3e} Im={row['lam_im']:.3f} "
                      f"L={row.get('L', '')} kappa={row.get('kappa', '')} P_outer={row['P_outer']:.1e} "
                      f"[{row['time']:.0f}s]", flush=True)
            bisect_thresholds(M, rows, pts, modes, m, n, shifts, kk, n_bisect)
            for r in rows:
                if "P_outer" not in r:
                    r["P_outer"] = outer_power_fraction(M, modes[round(r["b"], 6)][0]["w"], R)
            io_utils.write_csv(os.path.join(out, f"{name}_spectrum.csv"), rows)
            windows, thresholds = windows_and_thresholds(rows)
            summary = dict(n=n, m=m, Theta=T, R=R, h=M.h, b_min=sols[0]["b"], b_max=sols[-1]["b"],
                           n_points=len(rows), stop_up=info["stop_up"], stop_down=info["stop_down"],
                           thresholds=";".join(f"{a:.4f}-{b:.4f}({t})" for a, b, t in thresholds),
                           windows=";".join(f"{w[0]['b']:.3f}-{w[-1]['b']:.3f}" for w in windows),
                           max_P_outer=max(r["P_outer"] for r in rows))
            for iw, win in enumerate(windows):
                top = max(win, key=lambda r: r["growth"])
                summary.update({f"window{iw}_max_growth": top["growth"], f"window{iw}_b_at_max": top["b"],
                                f"window{iw}_lam_im_at_max": top["lam_im"], f"window{iw}_L_at_max": top.get("L", ""),
                                f"window{iw}_kappa_at_max": top.get("kappa", "")})
                s, mode = modes[round(top["b"], 6)]
                if mode is not None:
                    plot_mode(M, s, mode, top, os.path.join(out, f"{name}_mode_b{top['b']:.3f}.png"), name)
            if m == 2 and args.U_mark > 0:
                sU, rowU, modeU = spectrum_at_power(M, sols, args.U_mark, m, n, shifts, kk, R)
                summary["U_mark"] = args.U_mark
                if rowU is None:
                    print(f"   U={args.U_mark:g}: outside the family or not converged", flush=True)
                else:
                    summary.update(b_at_U_mark=sU["b"], growth_at_U_mark=rowU["growth"],
                                   lam_im_at_U_mark=rowU["lam_im"], unstable_at_U_mark=rowU["unstable"],
                                   L_at_U_mark=rowU.get("L", ""), kappa_at_U_mark=rowU.get("kappa", ""),
                                   P_outer_at_U_mark=rowU["P_outer"])
                    print(f"   U={args.U_mark:g}: b={sU['b']:.4f} growth={rowU['growth']:.3e} "
                          f"Im={rowU['lam_im']:.3f} L={rowU.get('L', '')} P_outer={rowU['P_outer']:.1e}", flush=True)
                    if rowU["unstable"] and modeU is not None:
                        plot_mode(M, sU, modeU, rowU, os.path.join(out, f"{name}_mode_U{args.U_mark:g}.png"), name)
            growth_plot(rows, os.path.join(out, f"{name}_growth.png"),
                        f"{name}: thresholds {summary['thresholds'] or 'none'}")
            with open(os.path.join(out, f"{name}_row.json"), "w", encoding="utf-8") as fh:
                json.dump(summary, fh, indent=1, default=io_utils.json_default)
            print(f"   thresholds: {summary['thresholds'] or 'none'}; windows: {summary['windows'] or 'none'}",
                  flush=True)
    if not args.job:
        assemble(out)
    print("\nE5c done.")


if __name__ == "__main__":
    main()
