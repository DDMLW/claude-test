"""
E3 —— 嬗变涡旋孤子的稳定性: 线性稳定性谱 + 加噪声传播
=====================================================
前提: 先用相同 --tag 跑完 E2a (本脚本读取 data/E2_transmutation_atlas/<tag>/ 中保存的解族)。

物理问题: 哪些 (n, m) 的嬗变涡旋是稳定的? 已知:
    圆柱形热样品中 m<=2 的环形涡旋稳定、m>=3 不稳定 (Kartashov et al., Opt. Express 15, 9378 (2007));
    矩形样品中 m=2 双涡旋稳定、m>=3 不稳定 (Ye et al., OL 35, 628 (2010))。
    多边形边界导致的嬗变结构 (中心反号奇点 + 卫星环) 能否让高荷涡旋稳定? —— 这是本课题的关键发现点,
    结果未知, 无论稳定与否都要如实记录。

方法:
  1) 线性稳定性: 全样品稀疏 shift-invert Arnoldi, 沿虚轴附近扫描 shift sigma = 0.02 + i omega
     (omega = 0 ... omega_max), 汇总所有找到的本征值, 增长率 = max Re(lambda)。
  2) 加噪声传播 (输入噪声方差 0.01, 与导师论文相同): 每个解族选代表点, 传播到 xi_max,
     记录偏离度、快照、奇点结构的演化。线性谱判不稳定、传播却长期保持 = 弱不稳定非线性饱和, 要单独标注。

输出 (results/E3_stability/<tag>/):
    <sample>_m<m>_spectrum.csv     每个 b: max Re(lambda), 最不稳定本征值
    <sample>_m<m>_stability.png    U(b) 曲线, 稳定段实线/不稳定段红叉, 及增长率 vs b
    <sample>_m<m>_prop_b<val>_*.png/csv  传播结果
    <sample>_m<m>_row.json         该解族的汇总行 (拆分并行时各任务分别写, 互不覆盖)
    stability_table.csv            总表 (由目录中所有 *_row.json 汇总; 拆分运行后用 --table_only 重建)

运行:
  python e3_stability.py --quick                       (约 10-20 分钟, 需先跑 E2a --quick)
  python e3_stability.py                               (数小时; 可按样品/荷拆开并行)
  python e3_stability.py --samples poly3 --ms 2 3      (只算一部分)
  python e3_stability.py --samples poly3 --ms 2 --bs 0.5 1.0 3.0 --tag check   (只在指定的 b 处算谱, 用于复核个别点)
  python e3_stability.py --table_only                  (拆分运行全部结束后重建 stability_table.csv)
"""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import build_sample, finalize_args, pick, sample_h, standard_args  # noqa: E402

from nlmodes import io_utils, plotting, propagate, stability  # noqa: E402
from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SAMPLES = ["rect", "poly3", "poly4", "poly5", "poly6", "disk"]
UNSTABLE_TOL = 1e-4      # max Re(lambda) 大于此值判为线性不稳定


def load_family(data_dir, kind, m, model):
    files = sorted(glob.glob(os.path.join(data_dir, f"{kind}_m{m}_b*.npz")))
    sols = []
    for f in files:
        d = io_utils.load_solution(f)
        w = io_utils.solution_on_model(d, model)
        sols.append(dict(w=w, b=d["b"], power=d["U"], n_t=None, converged=True))
    return sorted(sols, key=lambda s: s["b"])


def pick_spectrum_points(sols, stride):
    """
    选出要算谱的解: b 取在固定网格 b = k * stride * db 上 (db 为解族步长), 而不是简单地每隔 stride 个取一个,
    这样不同网格 (hq 与 hq_coarse) 的谱在相同的 b 处计算, 可以逐点比较。
    """
    if len(sols) < 3 or stride <= 1:
        return sols[::max(stride, 1)]
    bs = np.array([s["b"] for s in sols])
    db = float(np.median(np.diff(bs)))
    step = stride * db
    picked = [s for s in sols if abs(s["b"] / step - round(s["b"] / step)) < 0.2 / stride]
    return picked if len(picked) >= 2 else sols[::stride]


def assemble_table(out):
    """把目录中所有 <sample>_m<m>_row.json 汇总成 stability_table.csv (按样品、m 排序)。"""
    rows = []
    for f in glob.glob(os.path.join(out, "*_row.json")):
        with open(f, encoding="utf-8") as fh:
            rows.append(json.load(fh))
    order = {k: i for i, k in enumerate(SAMPLES)}
    rows.sort(key=lambda r: (order.get(r["sample"], 99), r["m"]))
    io_utils.write_csv(os.path.join(out, "stability_table.csv"), rows)
    print(f"stability_table.csv: {len(rows)} families")


def main():
    p = standard_args(__doc__)
    p.add_argument("--samples", nargs="+", default=SAMPLES)
    p.add_argument("--ms", nargs="+", type=int, default=None)
    p.add_argument("--rho", type=float, default=8.0)
    p.add_argument("--stride", type=int, default=None, help="每隔几个 b 算一次谱 (默认 2, quick 3)")
    p.add_argument("--bs", nargs="+", type=float, default=None,
                   help="只在这些 b 处算谱 (取最接近的已存解; 给出时忽略 --stride)")
    p.add_argument("--omega_max", type=float, default=None, help="shift 扫描的最大虚部 (默认 3, quick 2, hq 4)")
    p.add_argument("--zmax", type=float, default=None, help="传播距离 (默认 500, quick 40, hq 1000)")
    p.add_argument("--no_prop", action="store_true", help="只算线性谱, 不做传播")
    p.add_argument("--table_only", action="store_true", help="不计算, 只由已有的 *_row.json 重建总表")
    args = finalize_args(p, "E3_stability")
    out = args.res_dir
    if args.table_only:
        assemble_table(out)
        return
    e2_data = os.path.join(io_utils.PROJECT_ROOT, "data", "E2_transmutation_atlas", args.tag)
    if not os.path.isdir(e2_data):
        raise SystemExit(f"找不到 E2a 的数据目录 {e2_data}; 请先运行 E2a (相同的 --tag)")
    ms = args.ms or ([1, 2, 3, 4] if args.quick else [1, 2, 3, 4, 5, 6])
    stride = args.stride or (3 if args.quick else 2)
    om_max = args.omega_max or pick(args, quick=2.0, main=3.0, hq=4.0)
    d_om = pick(args, quick=0.25, main=0.15, hq=0.1)
    shifts = [0.02 + 1j * om for om in np.arange(0.0, om_max + 1e-9, d_om)]
    zmax = args.zmax or pick(args, quick=40.0, main=500.0, hq=1000.0)
    kk = pick(args, quick=12, main=16, hq=24)
    for kind in args.samples:
        size = (20.0, 40.0) if kind == "rect" else args.rho
        sample, nsym = build_sample(kind, sample_h(args, kind), size)
        M = NonlocalModel(sample)
        for m in ms:
            sols = load_family(e2_data, kind, m, M)
            if not sols:
                continue
            name = f"{kind}_m{m}"
            print(f"\n{name}: {len(sols)} stored solutions; "
                  + (f"spectra at b = {args.bs}" if args.bs else f"spectra every {stride}"))
            rows = []
            points = (sorted({id(s): s for s in (min(sols, key=lambda x: abs(x["b"] - bt)) for bt in args.bs)}
                                 .values(), key=lambda x: x["b"]) if args.bs else pick_spectrum_points(sols, stride))
            for s in points:
                st = stability.stability_arnoldi(M, s, shifts=shifts, k=kk)
                lam = st["eigenvalues"]
                j = int(np.argmax(lam.real))
                rows.append(dict(b=s["b"], U=s["power"], max_growth=st["max_growth"],
                                 lam_re=float(lam[j].real), lam_im=float(abs(lam[j].imag)),
                                 unstable=st["max_growth"] > UNSTABLE_TOL, n_eigs=len(lam),
                                 n_failed_shifts=st.get("n_failed_shifts", 0)))
                print(f"   b={s['b']:.3f} U={s['power']:.3f}  max Re(lambda) = {st['max_growth']:.3e} "
                      f"(Im = {abs(lam[j].imag):.3f})  [{st['time']:.0f}s]"
                      + (f"  WARNING: ARPACK not fully converged for {st['n_failed_shifts']} shifts"
                         if st.get("n_failed_shifts") else ""))
            io_utils.write_csv(os.path.join(out, f"{name}_spectrum.csv"), rows)
            fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.4))
            axs[0].plot([s["b"] for s in sols], [s["power"] for s in sols], "-", color="0.7")
            st_ = [r for r in rows if not r["unstable"]]
            un_ = [r for r in rows if r["unstable"]]
            axs[0].plot([r["b"] for r in st_], [r["U"] for r in st_], "o", color="C0", ms=4, label="stable")
            axs[0].plot([r["b"] for r in un_], [r["U"] for r in un_], "x", color="C3", ms=6, label="unstable")
            axs[0].set_xlabel("b")
            axs[0].set_ylabel("U")
            axs[0].legend(fontsize=7)
            axs[0].set_title(name)
            axs[1].semilogy([r["b"] for r in rows], [max(r["max_growth"], 1e-12) for r in rows], "o-")
            axs[1].axhline(UNSTABLE_TOL, color="k", ls=":")
            axs[1].set_xlabel("b")
            axs[1].set_ylabel("max Re(lambda)")
            fig.tight_layout()
            fig.savefig(os.path.join(out, f"{name}_stability.png"), dpi=130)
            plt.close(fig)

            # ---- 代表点传播: 最不稳定的点 与 最稳定的点 (若存在)
            prop = []
            if not args.no_prop and rows:
                picks = {min(rows, key=lambda r: r["max_growth"])["b"], max(rows, key=lambda r: r["max_growth"])["b"]}
                for bp in sorted(picks):
                    s = min(sols, key=lambda x: abs(x["b"] - bp))
                    snap = [0, zmax / 4, zmax / 2, 3 * zmax / 4, zmax]
                    res = propagate.propagate(M, propagate.add_noise(s["w"], 0.01, seed=7), zmax, dz=args.dz,
                                              diag_every=25, snapshot_zs=snap, w_ref=s["w"], nsym=nsym,
                                              verbose=False)
                    tag = f"{name}_prop_b{s['b']:.2f}"
                    ring = dg.ring_radius(dg.field_grid(M, s["w"]), M.grid, (0, 0), rmax=8.0)
                    plotting.plot_snapshots(M, res["snapshots"], os.path.join(out, tag + "_snapshots.png"),
                                            zoom=2 * ring + 2, title=f"{name}, b={s['b']:.2f}, noise 0.01",
                                            show_phase=True)
                    plotting.plot_timeseries(res, os.path.join(out, tag + "_timeseries.png"), title=tag)
                    dev = res["deviation"]
                    second_half = dev[len(dev) // 2:]
                    prop.append(dict(b=s["b"], zmax=zmax, max_dev=float(dev.max()),
                                     max_dev_second_half=float(second_half.max()),
                                     survived=bool(second_half.max() < 0.5)))
                    print(f"   propagation b={s['b']:.2f}: max deviation {dev.max():.3f} "
                          f"(second half {second_half.max():.3f})")
                io_utils.write_csv(os.path.join(out, f"{name}_propagation.csv"), prop)
            row = dict(sample=kind, n=nsym or 0, m=m, h=M.h, b_min=sols[0]["b"], b_max=sols[-1]["b"],
                       n_spectra=len(rows), n_unstable=sum(r["unstable"] for r in rows),
                       max_growth_overall=max(r["max_growth"] for r in rows),
                       stable_b=";".join(f"{r['b']:.2f}" for r in rows if not r["unstable"]),
                       propagation=";".join(f"b={x['b']:.2f}:{'ok' if x['survived'] else 'broke'}"
                                            for x in prop))
            with open(os.path.join(out, f"{name}_row.json"), "w", encoding="utf-8") as fh:
                json.dump(row, fh, indent=1, default=io_utils.json_default)
    if args.job:
        print("\nsplit job (--job given): run --table_only after all jobs finish to build stability_table.csv")
    else:
        assemble_table(out)
    print("\nE3 done.")


if __name__ == "__main__":
    main()
