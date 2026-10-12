"""
E3b —— 低功率失稳窗口: 阈值、失稳模式的结构、按角向通道分解的 Bogoliubov 谱
=========================================================================
背景 (Progress.md C8): hq 全扫描发现 m=2 涡旋在所有样品 (含圆盘) 中都有低功率振荡失稳窗口
(圆盘 b_c≈0.53), 三角形把窗口拉宽到 b_c≈2.7, 三角形中 m=1 也失稳。云端探索 (h=0.25) 给出的模式结构:
    圆盘 m=2      : L = -2 (u ~ e^{0 i theta}, v ~ e^{-4 i theta}) —— 双荷涡旋分裂为两个单位涡旋的四极模
    三角形 m=2    : b≈1.5 时 L = -1 与 +2 两个通道混合 (相差 n=3), 属 C_3 破缺类 kappa = -1;
                    低 b 时 L = -2 与 +1 混合 (kappa = +1)
    三角形 m=1    : L = -2 与 +1 混合 (kappa = +1)
假说: 圆盘中不同角向指数 L 的 Bogoliubov 通道互不耦合; C_n 边界把 L 与 L ± n 耦合, 原本不相干的、Krein 符号相反的
两支频率相遇时就会产生新的失稳窗口 (n=3 耦合最强, n=5,6 只有高阶耦合 -> 窗口与圆盘相同)。

本脚本对每个“情形” (样品:拓扑荷, 例如 poly3:2) 做:
  1) 细 b 网格解族 (b 步长 quick 0.2 / main 0.1 / hq 0.05);
  2) 每隔 0.1 (quick 0.2) 计算一次稳定性谱 (低功率 shift 扫描: omega = 0..1.2 步长 0.05, 外加 1.6..3.0 的稀疏 shift);
     对不稳定的点求最不稳定的 Bogoliubov 模并做角向谐波分析 (L、kappa、主导谐波);
  3) 稳定/不稳定交界处二分加密 (main 3 次, hq 4 次 -> b_c 精度约 0.006);
  4) 每个失稳窗口: 最不稳定点与最靠近阈值的点画模式图 (|u|、arg u、|v|、扰动后的光强);
  5) 最不稳定点沿本征模加小扰动后传播 (非线性结局: 分裂? 振荡? 破坏?)。
另有 --channels 模式 (只用于圆盘, --omega w1 w2 改变 shift 的频率区间): 每个 b 求 Im(lambda) <= 1 范围内的全部本征值并按角向通道 L 分类、给出 Krein 符号,
得到“圆盘的通道分解 Bogoliubov 谱”, 用来预言 C_n 样品中哪些通道 (L 与 L+n) 会相遇。

输出 (results/E3_stability/<tag>/):
    e3b_<样品>_m<m>_spectrum.csv     每个 b: 增长率、最不稳定本征值、L、kappa、主导谐波 (含二分加密的点)
    e3b_<样品>_m<m>_row.json         汇总行 (阈值区间、各窗口的模式类型) -> e3b_summary.csv (--table_only 重建)
    e3b_<样品>_m<m>_growth.png       增长率 vs b (标出 L/kappa)
    e3b_<样品>_m<m>_mode_b<val>.png  模式图;  e3b_<样品>_m<m>_evolution_*.png/csv  非线性结局
    e3b_disk_m<m>_channels.csv/png   (--channels) 圆盘的通道分解谱

运行:
  python e3b_instability_windows.py --quick                         (约 10 分钟: 圆盘与三角形 m=2, 粗网格)
  python e3b_instability_windows.py --cases poly3:2 --job poly3_m2  (一个情形; hq 单核约 2-4 h)
  python e3b_instability_windows.py --cases disk:2 --channels --job disk_m2_channels
  python e3b_instability_windows.py --table_only                    (拆分运行结束后重建 e3b_summary.csv)
"""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (build_sample, finalize_args, pick, rect_state_via_square, robust_family,  # noqa: E402
                    sample_h, standard_args, transfer_solution)
from stability_tools import (TOL, bisect_thresholds, low_power_shifts, plot_mode, spectrum_point,  # noqa: E402
                             windows_and_thresholds)

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, propagate, stability  # noqa: E402
from nlmodes.geometry import polygon_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ALL_CASES = ["rect:2", "poly3:2", "poly4:2", "poly5:2", "poly6:2", "disk:2", "poly3:1", "disk:1"]
B_HI = {("poly3", 2): 5.0, ("poly4", 2): 1.6}          # 其余情形默认 1.4 (三角形 m=2 到 5.0: 高 b 端的弱失稳, 见 README)


def parse_case(c):
    kind, m = c.split(":")
    return kind, int(m)


def build_family(M, kind, m, b_lo, b_hi, db, h):
    """
    扫描区间 [b_lo, b_hi] 内的解族。依次尝试:
      1) 在区间内起步 (b_start = min(3 + m/4, b_hi), 与 E2a 相同);
      2) 区间上端太低时 LG 初值常不收敛 (例如 20x40 矩形 m=2 在 b=1.4): 改在 b = 3 + m/4 起步再向下延拓, 只保留 b <= b_hi;
      3) 矩形: 与 E2a 相同, 先在正方形中求 b=5 的态, 沿宽度收缩到 20x40, 再沿 b 向下延拓。
    返回 (按 b 升序的解列表, info)。
    """
    shape = lambda W: guesses.lg_vortex(M, m, W)              # noqa: E731
    b_ref = 3.0 + 0.25 * m
    sols, info = robust_family(M, shape, b_lo, b_hi, db=db, b_start=min(b_ref, b_hi), verbose=False, align_db=True)
    if sols:
        return sols, info
    if b_ref > b_hi:
        sols, info = robust_family(M, shape, b_lo, b_ref, db=db, b_start=b_ref, verbose=False, align_db=True)
        sols = [s for s in sols if s["b"] <= b_hi + 1e-9]
        if sols:
            info.update(note=f"family started at b={b_ref:g} above the scan range")
            return sols, info
    if kind == "rect":
        Mr, s5 = rect_state_via_square(h, m, 5.0, verbose=False)
        if s5 is not None:
            w5 = transfer_solution(Mr, s5["w"], M)
            dn, info_dn = newton.continuation(M, w5, 5.0, b_lo, -db, verbose=False)
            sols = sorted([s for s in dn if s["b"] <= b_hi + 1e-9], key=lambda x: x["b"])
            if sols:
                return sols, dict(stop_up="", stop_down=info_dn["stop_reason"], b_min=sols[0]["b"],
                                  b_max=sols[-1]["b"], note="aspect-ratio continuation from the square at b=5")
    return [], info


def evolve_unstable(M, s, mode, row, nsym, args, out, name):
    """沿最不稳定模加小扰动后传播, 记录非线性结局。"""
    gam = max(row["growth"], 1e-3)
    zmax = float(np.clip(20.0 / gam, 50.0, pick(args, quick=60.0, main=300.0, hq=600.0)))
    pert = mode["u"] + np.conj(mode["v"])
    q0 = s["w"] + 1e-3 * np.abs(s["w"]).max() / np.abs(pert).max() * pert
    snap = list(np.linspace(0, zmax, 5))
    res = propagate.propagate(M, q0, zmax, dz=args.dz, diag_every=50, snapshot_zs=snap, w_ref=s["w"],
                              nsym=nsym or 4, track_singularities=True, sing_every=4, verbose=False)
    tag = f"{name}_evolution_b{s['b']:.2f}"
    ring = dg.ring_radius(dg.field_grid(M, s["w"]), M.grid, (0, 0))
    plotting.plot_snapshots(M, res["snapshots"], os.path.join(out, tag + "_snapshots.png"),
                            zoom=2.2 * ring + 1.0, title=f"{name} b={s['b']:.2f}: seeded along the unstable mode",
                            show_phase=True)
    plotting.plot_timeseries(res, os.path.join(out, tag + "_timeseries.png"), title=tag)
    io_utils.write_csv(os.path.join(out, tag + ".csv"),
                       [dict(z=z, deviation=d, central=c, n_singularities=(len(sg) if sg is not None else ""))
                        for z, d, c, sg in zip(res["z"], res["deviation"], res["central"], res["singularities"])])
    dev = res["deviation"]
    return dict(evolution_zmax=zmax, evolution_max_dev=float(dev.max()),
                evolution_final_central=float(res["central"][-1]),
                evolution_saturated=bool(dev[len(dev) // 2:].max() < 2 * dev[len(dev) // 2]))


def channels(M, sols, m, args, out, name, stride):
    """圆盘: Im(lambda) <= 1.0 (或 --omega 给出的区间) 内的全部本征值, 按角向通道 L 分类并给出 Krein 符号。"""
    w1, w2 = args.omega or (0.05, 1.0)
    sigmas = [0.0 + 1j * w for w in np.arange(w1, w2 + 0.01, pick(args, quick=0.25, main=0.15, hq=0.1))]
    kk = pick(args, quick=8, main=12, hq=16)
    rows = []
    for s in sols[::stride]:
        seen = []
        for sg in sigmas:
            for lam, md in stability.eigenpairs_near(M, s, sg, k=kk):
                if lam.imag < 0.01 or any(abs(lam - x) < 1e-6 for x in seen):
                    continue
                seen.append(lam)
                ac = stability.mode_angular_content(M, md, m, w=s["w"])
                rows.append(dict(b=s["b"], lam_re=lam.real, lam_im=lam.imag, L=ac["L"], L_from_v=ac["L_from_v"],
                                 krein=stability.krein_sign(md), w_u=round(ac["w_u"], 3)))
        print(f"   channels b={s['b']:.2f}: {len(seen)} eigenvalues", flush=True)
    io_utils.write_csv(os.path.join(out, f"{name}_channels.csv"), rows)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    # 固定配色: L = -3..3 各一种颜色 (色盲安全的固定顺序), |L| > 3 归入灰色 "other";
    # Krein 符号用形状区分 (实心圆 +, 空心三角 -), 颜色不是唯一线索
    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
    for L in list(range(-3, 4)) + [None]:
        col = palette[L + 3] if L is not None else "0.6"
        sel_L = [r for r in rows if (r["L"] == L if L is not None else abs(r["L"]) > 3)]
        for kr, mk, mfc in ((1, "o", col), (-1, "^", "none")):
            sel = [r for r in sel_L if r["krein"] == kr]
            if sel:
                ax.plot([r["b"] for r in sel], [r["lam_im"] for r in sel], mk, ms=4.5, ls="none", color=col,
                        mfc=mfc, label=(f"L={L:+d}" if L is not None else "|L|>3") if kr == 1 or not
                        [r for r in sel_L if r["krein"] == 1] else None)
    uns = [r for r in rows if r["lam_re"] > TOL]
    ax.plot([r["b"] for r in uns], [r["lam_im"] for r in uns], "s", ms=9, mfc="none", mec="k", ls="none",
            label="unstable")
    ax.set_xlabel("b")
    ax.set_ylabel("Im(lambda)")
    ax.set_title(f"{name}: Bogoliubov spectrum by angular channel L (o: Krein +, ^: Krein -)", fontsize=8)
    ax.legend(fontsize=6, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{name}_channels.png"), dpi=130)
    plt.close(fig)


def assemble(out):
    rows = []
    for f in glob.glob(os.path.join(out, "e3b_*_row.json")):
        with open(f, encoding="utf-8") as fh:
            rows.append(json.load(fh))
    rows.sort(key=lambda r: (r["m"], r["sample"], r.get("theta0_deg", 0.0), r["b_min"]))
    io_utils.write_csv(os.path.join(out, "e3b_summary.csv"), rows)
    print(f"e3b_summary.csv: {len(rows)} cases")


def main():
    p = standard_args(__doc__)
    p.add_argument("--cases", nargs="+", default=None, help="样品:拓扑荷 列表, 默认全部 (quick: disk:2 poly3:2)")
    p.add_argument("--rho", type=float, default=8.0)
    p.add_argument("--b_range", nargs=2, type=float, default=None, help="只算这个 b 区间 (例如阈值附近的加密复核)")
    p.add_argument("--theta0", type=float, default=0.0,
                   help="多边形相对网格的旋转角 (度); 用于检验阈值是否受网格 C_4 各向异性影响 (旋转检验)")
    p.add_argument("--channels", action="store_true", help="只做圆盘的通道分解谱 (cases 必须是 disk:m)")
    p.add_argument("--omega", nargs=2, type=float, default=None,
                   help="--channels 模式: 只在 Im(lambda) 的这个区间内放 shift (默认 0.05-1.0; 高 b 或大半径时频率更高/更低)")
    p.add_argument("--no_evolution", action="store_true", help="不做非线性结局的传播")
    p.add_argument("--spec_step", type=float, default=None,
                   help="每隔多少 b 算一次谱 (默认 quick 0.2 / main 0.1 / hq 0.1; 细网格加密任务用 0.2 省时, 阈值靠二分)")
    p.add_argument("--table_only", action="store_true", help="不计算, 只由 e3b_*_row.json 重建 e3b_summary.csv")
    args = finalize_args(p, "E3_stability")
    out = args.res_dir
    if args.table_only:
        assemble(out)
        return
    cases = [parse_case(c) for c in (args.cases or (["disk:2", "poly3:2"] if args.quick else ALL_CASES))]
    db = pick(args, quick=0.2, main=0.1, hq=0.05)
    spec_step = args.spec_step or pick(args, quick=0.2, main=0.1, hq=0.1)
    stride = max(1, int(round(spec_step / db)))
    shifts = low_power_shifts(args)
    kk = pick(args, quick=10, main=14, hq=18)
    n_bisect = pick(args, quick=1, main=3, hq=4)
    for kind, m in cases:
        name = f"e3b_{kind}_m{m}" + (f"_rot{args.theta0:g}" if args.theta0 else "")
        size = (20.0, 40.0) if kind == "rect" else args.rho
        if args.theta0 and kind.startswith("poly"):
            sample, nsym = polygon_sample(int(kind[4:]), sample_h(args, kind), rho=size,
                                          theta0=np.radians(args.theta0)), int(kind[4:])
        else:
            sample, nsym = build_sample(kind, sample_h(args, kind), size)
        M = NonlocalModel(sample)
        b_lo, b_hi = args.b_range or (0.05 if not args.quick else 0.2, B_HI.get((kind, m), 1.4))
        print(f"\n{M.describe()}\n{name}: b in [{b_lo}, {b_hi}], db={db}, spectra every {spec_step}", flush=True)
        sols, info = build_family(M, kind, m, b_lo, b_hi, db, sample_h(args, kind))
        print(f"   family: {len(sols)} solutions, b in [{sols[0]['b'] if sols else None}, "
              f"{sols[-1]['b'] if sols else None}] {info.get('note', '')}", flush=True)
        if not sols:
            continue
        if args.channels:
            channels(M, sols, m, args, out, f"{kind}_m{m}" if kind != "disk" else f"e3b_disk_m{m}", stride)
            continue
        pts = sols[::stride]
        rows, modes = [], {}
        for s in pts:
            row, mode = spectrum_point(M, s, m, nsym, shifts, kk)
            rows.append(row)
            modes[round(s["b"], 6)] = (s, mode)
            print(f"   b={s['b']:.3f} growth={row['growth']:.3e} Im={row['lam_im']:.3f} L={row.get('L', '')} "
                  f"kappa={row.get('kappa', '')} [{row['time']:.0f}s]", flush=True)
        bisect_thresholds(M, rows, pts, modes, m, nsym, shifts, kk, n_bisect)
        io_utils.write_csv(os.path.join(out, f"{name}_spectrum.csv"), rows)
        windows, thresholds = windows_and_thresholds(rows)
        summary = dict(sample=kind, m=m, n=nsym or 0, h=M.h, theta0_deg=args.theta0, b_min=sols[0]["b"],
                       b_max=sols[-1]["b"], n_points=len(rows),
                       thresholds=";".join(f"{a:.4f}-{b:.4f}({t})" for a, b, t in thresholds),
                       windows=";".join(f"{w[0]['b']:.3f}-{w[-1]['b']:.3f}" for w in windows))
        # ---- 每个窗口: 模式图 + 非线性结局
        for iw, win in enumerate(windows):
            top = max(win, key=lambda r: r["growth"])
            edge = win[-1] if win[-1]["b"] > top["b"] else win[0]
            summary[f"window{iw}_max_growth"] = top["growth"]
            summary[f"window{iw}_b_at_max"] = top["b"]
            summary[f"window{iw}_L_at_max"] = top.get("L", "")
            summary[f"window{iw}_kappa_at_max"] = top.get("kappa", "")
            summary[f"window{iw}_L_at_edge"] = edge.get("L", "")
            for r in {id(top): top, id(edge): edge}.values():
                s, mode = modes[round(r["b"], 6)]
                if mode is not None:
                    plot_mode(M, s, mode, r, os.path.join(out, f"{name}_mode_b{r['b']:.3f}.png"), name)
            if not args.no_evolution:
                s, mode = modes[round(top["b"], 6)]
                if mode is not None:
                    summary.update({f"window{iw}_{k}": v for k, v in
                                    evolve_unstable(M, s, mode, top, nsym, args, out, name).items()})
        # ---- 增长率图
        fig, ax = plt.subplots(figsize=(6, 3.6))
        ax.plot([r["b"] for r in rows], [max(r["growth"], 0.0) for r in rows], "o-", ms=3, color="#2a78d6")
        for r in rows:
            if r["unstable"] and r.get("L", "") != "":
                ax.annotate(f"{r['L']:+d}", (r["b"], r["growth"]), fontsize=6, xytext=(0, 4),
                            textcoords="offset points", ha="center")
        ax.axhline(TOL, color="0.5", lw=0.8)
        ax.set_xlabel("b")
        ax.set_ylabel("max Re(lambda)  (labels: angular index L)")
        ax.set_title(f"{name}: thresholds {summary['thresholds']}", fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(out, f"{name}_growth.png"), dpi=130)
        plt.close(fig)
        with open(os.path.join(out, f"{name}_row.json"), "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=1, default=io_utils.json_default)
        print(f"   thresholds: {summary['thresholds'] or 'none'}; windows: {summary['windows'] or 'none'}", flush=True)
    if not args.job:
        assemble(out)
    print("\nE3b done.")


if __name__ == "__main__":
    main()
