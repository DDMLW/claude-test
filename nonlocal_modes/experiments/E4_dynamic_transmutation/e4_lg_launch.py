"""
E4 —— 动态嬗变: 把实验上最容易产生的正则 LG 涡旋光束 (螺旋相位板) 射入多边形热样品
==================================================================================
这是“可被实验验证”的预言: 输入光束是标准的 m 阶涡旋 (单个中心奇点), 样品材料均匀,
只有截面形状是正 n 边形。理论预言:
    * 角赝动量 mu = m mod n 严格守恒 (C_n 不可约表示的权重在传播中不变)
    * 传播中中心奇点自发“嬗变”为 v (例如 n=3, m=2: 中心 +2 -> -1, 同时放出 3 个 +1 卫星)
    * 圆盘样品 (对照组) 中不发生嬗变 (仅有网格造成的极微小效应)

输出 (results/E4_dynamic_transmutation/<tag>/):
    <sample>_m<m>_U<U>_snapshots.png   振幅+相位快照 (标出奇点)
    <sample>_m<m>_U<U>_worldlines.png  奇点到中心距离 r 随 xi 的演化 (“世界线”)
    <sample>_m<m>_U<U>_irreps.png/csv  C_n 不可约表示权重随 xi (守恒检验)
    <sample>_m<m>_U<U>_interferogram.png  输出面与倾斜平面波的干涉图 (实验中看到的“叉形条纹”)
    e4_summary.csv                     每个 run 的最终中心荷 (小圆绕数)、卫星数、mu 守恒误差
                                       (由目录中所有 *_row.json 汇总; 拆分运行后可用 --table_only 重建)

运行:
  python e4_lg_launch.py --quick    (约 5-10 分钟)
  python e4_lg_launch.py            (约 1-3 小时)
  python e4_lg_launch.py --samples poly3 --job poly3     (只算一个样品, 便于并行)
  python e4_lg_launch.py --table_only                    (拆分运行全部结束后重建 e4_summary.csv)
"""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import build_sample, finalize_args, pick, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, plotting, propagate, theory  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def interferogram(model, q, kx=4.0, ref_amp=None):
    """与倾斜平面波 A exp(i kx x) 的干涉光强 (实验中用来读出相位奇点: 每个 +1 奇点对应一个叉)。"""
    qg = dg.field_grid(model, q)
    A = ref_amp if ref_amp is not None else 0.5 * np.abs(qg).max()
    return np.abs(qg + A * np.exp(1j * kx * model.grid.X) * (np.abs(qg) > 1e-3 * np.abs(qg).max())) ** 2


def assemble_summary(out):
    """把目录中所有 <sample>_m<m>_U<U>_row.json 汇总成 e4_summary.csv。"""
    rows = []
    for f in glob.glob(os.path.join(out, "*_row.json")):
        with open(f, encoding="utf-8") as fh:
            rows.append(json.load(fh))
    rows.sort(key=lambda r: (r["sample"], r["m"], r["U"]))
    io_utils.write_csv(os.path.join(out, "e4_summary.csv"), rows)
    print(f"e4_summary.csv: {len(rows)} runs")


def main():
    p = standard_args(__doc__)
    p.add_argument("--samples", nargs="+", default=["poly3", "poly4", "poly5", "poly6", "disk"])
    p.add_argument("--ms", nargs="+", type=int, default=None, help="默认 2 3 4 (quick: 2 3)")
    p.add_argument("--powers", nargs="+", type=float, default=None, help="输入能流 U 列表 (默认 10 20)")
    p.add_argument("--rho", type=float, default=8.0)
    p.add_argument("--w0", type=float, default=2.0, help="LG 光束宽度参数")
    p.add_argument("--zmax", type=float, default=None, help="传播距离 (默认 200, quick 30, hq 400)")
    p.add_argument("--noise", type=float, default=0.0, help="输入噪声方差 (默认 0, 可设 0.01 检验鲁棒性)")
    p.add_argument("--table_only", action="store_true", help="不计算, 只由已有的 *_row.json 重建 e4_summary.csv")
    args = finalize_args(p, "E4_dynamic_transmutation")
    out = args.res_dir
    if args.table_only:
        assemble_summary(out)
        return
    ms = args.ms or ([2, 3] if args.quick else [2, 3, 4])
    powers = args.powers or ([15.0] if args.quick else [10.0, 20.0])
    zmax = args.zmax or pick(args, quick=30.0, main=200.0, hq=400.0)
    nsnap = 6
    for kind in args.samples:
        sample, nsym = build_sample(kind, args.h, args.rho)
        M = NonlocalModel(sample)
        print("\n" + M.describe())
        for m in ms:
            for U in powers:
                q0 = guesses.lg_vortex(M, m, args.w0)
                q0 *= np.sqrt(U / M.power(q0))
                if args.noise > 0:
                    q0 = propagate.add_noise(q0, args.noise, seed=11)
                snap = list(np.linspace(0, zmax, nsnap))
                res = propagate.propagate(M, q0, zmax, dz=args.dz, diag_every=25, snapshot_zs=snap,
                                          nsym=nsym or 4, track_singularities=True, verbose=False)
                tag = f"{kind}_m{m}_U{U:g}"
                ring = dg.ring_radius(dg.field_grid(M, q0), M.grid, (0, 0), rmax=8.0)
                plotting.plot_snapshots(M, res["snapshots"], os.path.join(out, tag + "_snapshots.png"),
                                        zoom=2 * ring + 1.5, title=f"LG m={m} into {kind}, U={U:g}",
                                        show_phase=True)
                # 奇点世界线
                fig, ax = plt.subplots(figsize=(5, 3.4))
                for z, sg in zip(res["z"], res["singularities"]):
                    for (x, y, c) in sg:
                        r = np.hypot(x, y)
                        if r < 2 * ring:
                            ax.plot(z, r, ".", ms=2, color="red" if c > 0 else "blue")
                ax2 = ax.twinx()
                ax2.plot(res["z"], res["central"], "k-", lw=1)
                ax2.set_ylabel("central charge (black line)")
                ax.set_xlabel("xi")
                ax.set_ylabel("distance of unit singularity from centre")
                ax.set_title(f"{tag}: red +1, blue -1")
                fig.tight_layout()
                fig.savefig(os.path.join(out, tag + "_worldlines.png"), dpi=130)
                plt.close(fig)
                # 不可约表示权重
                irr = res["irreps"]
                mus = sorted(irr[0])
                io_utils.write_csv(os.path.join(out, tag + "_irreps.csv"),
                                   [dict(z=z, **{f"mu{k}": v[k] for k in mus}) for z, v in zip(res["z"], irr)])
                fig, ax = plt.subplots(figsize=(5, 3.2))
                for k in mus:
                    ax.semilogy(res["z"], [max(v[k], 1e-16) for v in irr], label=f"mu={k}")
                ax.set_xlabel("xi")
                ax.set_ylabel(f"C_{nsym or 4} irrep weight")
                ax.legend(fontsize=7)
                fig.tight_layout()
                fig.savefig(os.path.join(out, tag + "_irreps.png"), dpi=130)
                plt.close(fig)
                # 干涉图
                fig, axs = plt.subplots(1, 2, figsize=(7, 3.3))
                for ax, (z, key) in zip(axs, ((0.0, 0.0), (zmax, zmax))):
                    qz = res["snapshots"][min(res["snapshots"], key=lambda zz: abs(zz - key))]
                    ax.imshow(interferogram(M, qz), origin="lower", extent=M.grid.extent, cmap="gray")
                    ax.set_xlim(-2 * ring - 1, 2 * ring + 1)
                    ax.set_ylim(-2 * ring - 1, 2 * ring + 1)
                    ax.set_title(f"interferogram xi={z:g}")
                fig.tight_layout()
                fig.savefig(os.path.join(out, tag + "_interferogram.png"), dpi=130)
                plt.close(fig)
                # 末态结构
                qf = dg.field_grid(M, res["final"])
                Wc = dg.winding_on_circle(qf, M.grid, (0, 0), 0.5 * M.h)
                # 只统计亮环以内 (r < 0.8 r_ring) 的奇点 = 嬗变卫星; 外围辐射干涉产生的正负奇点对不计入
                ring_f = dg.ring_radius(qf, M.grid, (0, 0), rmax=8.0)
                reg = dg.beam_region(qf, M.grid, rel_threshold=1e-3) & M.opt_mask & (M.grid.R < 0.8 * ring_f)
                st = dg.singularity_structure(dg.find_singularities(qf, M.grid, reg, exclude_center=(0, 0),
                                                                    exclude_radius=1.25 * M.h),
                                              (0, 0), 1.25 * M.h)
                mu0 = theory.irrep_label(m, nsym or 4)
                drift = max(abs(v[mu0] - irr[0][mu0]) for v in irr)
                pred = theory.selection_rule(m, nsym)["central"] if nsym else m
                row = dict(sample=kind, n=nsym or 0, m=m, U=U, h=M.h, zmax=zmax, noise=args.noise,
                           central_pred=pred, central_final=round(Wc, 2),
                           n_satellites_final=st["n_satellites"], r_s_final=st["satellite_r_mean"],
                           mu=mu0, mu_weight_initial=irr[0][mu0], mu_weight_drift=drift,
                           power_drift=float(np.ptp(res["power"]) / res["power"][0]))
                with open(os.path.join(out, tag + "_row.json"), "w", encoding="utf-8") as fh:
                    json.dump(row, fh, indent=1, default=io_utils.json_default)
                print(f"  m={m} U={U:g}: central (pred {pred}) final {Wc:+.2f}, satellites {st['n_satellites']}, "
                      f"mu-weight drift {drift:.1e}  [{res['time']:.0f}s]")
    if args.job:
        print("\nsplit job (--job given): run --table_only after all jobs finish to build e4_summary.csv")
    else:
        assemble_summary(out)
    print("\nE4 done.")


if __name__ == "__main__":
    main()
