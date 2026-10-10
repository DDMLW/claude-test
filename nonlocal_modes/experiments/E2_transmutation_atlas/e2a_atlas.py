"""
E2a —— 拓扑荷嬗变“周期表” (transmutation atlas): 本课题的核心数值结果
=====================================================================
物理问题: 拓扑荷为 m 的涡旋孤子放在旋转对称阶数为 n 的热样品中 (材料完全均匀, 对称性只来自边界),
其相位奇点结构是什么? 理论 (nlmodes/theory.py, docs/theory.md) 预言:
    * 角赝动量 mu = m mod n 守恒, 解只含角谐波 l = m - j n
    * 中心奇点荷 v = |l| 最小的谐波 (例: n=3,m=2 -> v=-1; n=4,m=3 -> v=-1; n=2,m=2 -> v=0)
    * 其余 m - v 个单位荷以 C_n 对称的“卫星环”分布, 环的个数/半径标度由牛顿多边形给出
    * n | 2m 时 (+m 与 -m 被直接耦合), 涡旋只在功率阈值以上存在 (矩形中的双涡旋就是一例)

本实验: 对样品 {rect 20x40 (n=2), poly3, poly4, poly5, poly6 (保角半径 rho=8), disk (R=8, 对照)}
与 m = 1..6, 沿 b 求完整解族, 记录奇点结构并与理论逐项对比。

输出 (results/E2_transmutation_atlas/<tag>/):
    atlas_amplitude.png, atlas_phase.png   6x6 图谱 (行: 样品, 列: m), b = b_atlas 处
    atlas_table.csv                        理论 vs 数值 (中心荷、卫星数、mu、是否一致)
    <sample>_m<m>_summary.csv              每个解族每个 b 的全部诊断量
    <sample>_m<m>_b<val>.png               代表性模式图

运行:
  python e2a_atlas.py --quick                          (约 10 分钟: 粗网格, m<=4)
  python e2a_atlas.py                                  (约 1-2 小时)
  python e2a_atlas.py --samples poly3 poly4 --ms 2 3   (只算一部分, 便于多开几个终端并行)
  python e2a_atlas.py --plot_only                      (拆分运行全部结束后, 从 data/ 重建完整的表与图谱)
"""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (build_sample, finalize_args, pick_closest, rect_state_via_square,  # noqa: E402
                    robust_family, save_family, standard_args, transfer_solution)

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, theory  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SAMPLES = ["rect", "poly3", "poly4", "poly5", "poly6", "disk"]


def sample_size(kind, rho):
    return (20.0, 40.0) if kind == "rect" else rho


def compute_families(args, kind, ms, b_lo, b_hi, db):
    """计算一个样品的全部解族, 原始解写入 data/, 摘要写入 results/; 另写 <kind>_m<m>_info.json 记录端点信息。"""
    out, dat = args.res_dir, args.data_dir
    sample, nsym = build_sample(kind, args.h, sample_size(kind, args.rho))
    M = NonlocalModel(sample)
    print("\n" + M.describe())
    for m in ms:
        name = f"{kind}_m{m}"
        for old in glob.glob(os.path.join(dat, f"{name}_b*.npz")):
            os.remove(old)                      # 重新计算时清除旧数据, 避免新旧混杂
        sols, info = robust_family(M, lambda W, m=m: guesses.lg_vortex(M, m, W), b_lo, b_hi, db=db,
                                   b_start=min(3.0 + 0.25 * m, b_hi), verbose=False)
        if not sols and kind == "rect":
            print(f"  m={m}: direct solve failed; aspect-ratio continuation from the 40x40 square at b=5")
            Mr, s5 = rect_state_via_square(args.h, m, 5.0)
            if s5 is not None:
                w5 = transfer_solution(Mr, s5["w"], M)
                up, _ = newton.continuation(M, w5, 5.0, b_hi, db, verbose=False)
                dn, info_dn = newton.continuation(M, w5, 5.0, b_lo, -db, verbose=False)
                sols = sorted(dn[1:] + up, key=lambda x: x["b"])
                info = dict(stop_down=info_dn["stop_reason"], b_min=sols[0]["b"] if sols else None,
                            b_max=sols[-1]["b"] if sols else None)
        print(f"  m={m}: {len(sols)} solutions; b in [{info.get('b_min')}, {info.get('b_max')}]; "
              f"low end: {info['stop_down']}")
        with open(os.path.join(dat, f"{name}_info.json"), "w", encoding="utf-8") as f:
            json.dump(dict(n_solutions=len(sols), **{k: info.get(k) for k in ("b_min", "b_max", "stop_down")}),
                      f, ensure_ascii=False)
        if sols:
            save_family(M, sols, out, dat, name, nsym=nsym, m=m, plot_bs=[0.5, 1.0, 2.0, 4.0])


def load_family(dat, kind, m, model):
    sols = []
    for f in sorted(glob.glob(os.path.join(dat, f"{kind}_m{m}_b*.npz"))):
        d = io_utils.load_solution(f)
        sols.append(dict(w=io_utils.solution_on_model(d, model), b=d["b"], power=d["U"], n_t=None,
                         converged=True))
    return sorted(sols, key=lambda s: s["b"])


def assemble(args, ms):
    """从 data/ 中已保存的所有解族汇总出 atlas_table.csv 与图谱 (拆分运行后用 --plot_only 重建完整图谱)。"""
    out, dat = args.res_dir, args.data_dir
    table, atlas = [], {}
    for kind in SAMPLES:
        info_files = [os.path.join(dat, f"{kind}_m{m}_info.json") for m in ms]
        if not any(os.path.exists(f) for f in info_files):
            continue
        sample, nsym = build_sample(kind, args.h, sample_size(kind, args.rho))
        M = NonlocalModel(sample)
        for m, finfo in zip(ms, info_files):
            if not os.path.exists(finfo):
                continue
            with open(finfo, encoding="utf-8") as f:
                info = json.load(f)
            sols = load_family(dat, kind, m, M)
            if not sols:
                table.append(dict(sample=kind, n=nsym or 0, m=m, status="no solution found"))
                continue
            s_at = pick_closest(sols, "b", args.b_atlas)
            atlas[(kind, m)] = (M, s_at)
            r_at = dg.mode_summary(M, s_at, nsym=nsym, m=m)
            pred = theory.selection_rule(m, nsym) if nsym else dict(mu=m, central=m, n_satellites=0,
                                                                       degenerate=False, rings=[])
            num_central = int(round(r_at["W_small"]))
            resolved = (pred["central"] == m) or (r_at.get("r_s_fit", 0) >= 2 * args.h)
            agree = (num_central == pred["central"]) and abs(r_at["total_charge"] - m) < 0.1
            table.append(dict(
                sample=kind, n=nsym or 0, m=m, b=round(s_at["b"], 4), U=round(s_at["power"], 4),
                mu_pred=pred["mu"], mu_num=r_at.get("mu_dominant", m), mu_weight=r_at.get("mu_weight", 1.0),
                central_pred=pred["central"], central_num=num_central,
                satellites_pred=pred["n_satellites"], satellite_charge_num=int(r_at["satellite_charge"]),
                n_satellites_detected=r_at["n_satellites"], r_s_detected=r_at["satellite_r_mean"],
                r_s_fit=r_at.get("r_s_fit"), total_charge=r_at["total_charge"],
                W_profile=r_at["W_profile"], unresolved_charge=r_at["unresolved_charge"],
                resolved=resolved, degenerate_pred=pred["degenerate"], b_min_reached=info.get("b_min"),
                low_end=info.get("stop_down"),
                match=agree if resolved else ("unresolved (agrees)" if agree else "unresolved")))
            print(f"  {kind} m={m} b={s_at['b']:.2f}: central pred/num = {pred['central']}/{num_central}, "
                  f"satellites pred/num = {pred['n_satellites']}/{r_at['n_satellites']}, "
                  f"r_s(det/fit) = {r_at['satellite_r_mean']:.3f}/{r_at.get('r_s_fit', 0):.3f}")
    io_utils.write_csv(os.path.join(out, "atlas_table.csv"), table)

    kinds = [k for k in SAMPLES if any((k, m) in atlas for m in ms)]
    if atlas:
        for what in ("amplitude", "phase"):
            fig, axs = plt.subplots(len(kinds), len(ms), figsize=(2.3 * len(ms), 2.3 * len(kinds)),
                                    squeeze=False)
            for i, kind in enumerate(kinds):
                for j, m in enumerate(ms):
                    ax = axs[i, j]
                    ax.set_xticks([])
                    ax.set_yticks([])
                    if i == 0:
                        ax.set_title(f"m={m}", fontsize=9)
                    if j == 0:
                        ax.set_ylabel(kind, fontsize=9)
                    if (kind, m) not in atlas:
                        ax.text(0.5, 0.5, "no solution\nin b range", ha="center", va="center", fontsize=7,
                                transform=ax.transAxes)
                        continue
                    M, s = atlas[(kind, m)]
                    q = dg.field_grid(M, s["w"])
                    g = M.grid
                    zoom = 1.6 * dg.ring_radius(q, g, (0, 0), rmax=8.0) + 1.0
                    if what == "amplitude":
                        ax.imshow(np.abs(q), origin="lower", extent=g.extent, cmap="inferno")
                    else:
                        ph = np.where(np.abs(q) > 1e-4 * np.abs(q).max(), np.angle(q), np.nan)
                        ax.imshow(ph, origin="lower", extent=g.extent, cmap="twilight", vmin=-np.pi, vmax=np.pi)
                        sg, central = plotting.singularities_for_plot(q, M, (0.0, 0.0))
                        for (x, y, c) in sg:
                            ax.plot(x, y, "o", ms=3.5, mfc="none", mec="red" if c > 0 else "cyan", mew=1.2)
                        plotting.annotate_center(ax, (0.0, 0.0), central, fontsize=7)
                    ax.set_xlim(-zoom, zoom)
                    ax.set_ylim(-zoom, zoom)
            fig.suptitle(f"transmutation atlas ({what}), b ~ {args.b_atlas:g}")
            fig.tight_layout()
            fig.savefig(os.path.join(out, f"atlas_{what}.png"), dpi=140)
            plt.close(fig)
    n_res = [r for r in table if r.get("resolved")]
    n_match = sum(1 for r in n_res if r.get("match") is True)
    print(f"\ntheory/numerics agreement: {n_match}/{len(n_res)} resolved families; "
          f"{len(table) - len(n_res)} unresolved/no-solution (see atlas_table.csv)")


def main():
    p = standard_args(__doc__)
    p.add_argument("--samples", nargs="+", default=SAMPLES, help="要计算的样品 (默认全部)")
    p.add_argument("--ms", nargs="+", type=int, default=None, help="拓扑荷列表 (默认 1..6, quick 1..4)")
    p.add_argument("--rho", type=float, default=8.0, help="多边形保角半径 / 圆盘半径 (默认 8)")
    p.add_argument("--b_atlas", type=float, default=1.0, help="图谱所用的传播常数 (默认 1.0)")
    p.add_argument("--plot_only", action="store_true",
                   help="不计算, 只从 data/ 中已有的解族重建 atlas_table.csv 与图谱 (拆分运行后使用)")
    args = finalize_args(p, "E2_transmutation_atlas")
    ms = args.ms or ([1, 2, 3, 4] if args.quick else [1, 2, 3, 4, 5, 6])
    b_lo, b_hi, db = (0.5, 6.0, 0.5) if args.quick else (0.2, 7.0, 0.2)
    if not args.plot_only:
        for kind in args.samples:
            compute_families(args, kind, ms, b_lo, b_hi, db)
    print("\nassembling atlas from all families stored in", args.data_dir)
    assemble(args, ms)
    print("E2a done.")


if __name__ == "__main__":
    main()
