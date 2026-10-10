"""
E2b —— 卫星奇点环半径的标度律 (定量检验“牛顿多边形”理论)
=========================================================
理论: 正 n 边形中, 中心处光束感受的各向异性强度 eps ∝ U / rho^n (E1a 已验证系数)。
牛顿多边形给出卫星环半径 r_s ∝ eps^s (s 见 theory.selection_rule 的 rings), 因此在固定能流 U 下

        r_s ∝ rho^(-n s)        (对数-对数斜率 = -n s)

    (n,m) = (3,2): v=-1, s=1    -> 斜率 -3
    (n,m) = (4,3): v=-1, s=1/2  -> 斜率 -2
    (n,m) = (3,3): v= 0, s=1/3  -> 斜率 -1
    (n,m) = (4,4): v= 0, s=1/4  -> 斜率 -1
    (n=5,6 的物理各向异性在多边形中太弱、会被网格误差掩盖, 其标度律改在 E5a 中用加热图样检验)
    (n,m) = (2,2) 矩形 (双涡旋): s=1/2, 但矩形没有单一尺寸参数, 不在此测

方法: 定功率牛顿迭代 + 沿 rho 的参数延拓。从大样品 (弱各向异性, 正则 LG 初值几乎就是解) 出发,
逐步缩小 rho; 每个 rho 记录 r_s (网格检测值 与 中心谐波拟合值 r_s_fit, 后者可小于网格间距)。

输出: scaling_<n>_<m>.csv, scaling.png (log-log, 带理论斜率参考线), scaling_fits.csv (拟合斜率 vs 理论)

运行:
  python e2b_scaling.py --quick     (约 5 分钟)
  python e2b_scaling.py             (约 1-2 小时; 建议 --h 0.15 做一次加密复核)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, initial_solution, param_continuation, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, theory  # noqa: E402
from nlmodes.geometry import polygon_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CASES = [(3, 2, 12.0), (4, 3, 18.0), (3, 3, 18.0), (4, 4, 24.0)]     # (n, m, 固定能流 U)


def main():
    p = standard_args(__doc__)
    p.add_argument("--cases", nargs="+", default=None, help="形如 3,2 4,3 的 (n,m) 列表")
    args = finalize_args(p, "E2_transmutation_atlas")
    out = args.res_dir
    cases = CASES
    if args.cases:
        want = {tuple(int(v) for v in c.split(",")) for c in args.cases}
        cases = [c for c in CASES if (c[0], c[1]) in want]
    rhos = [16.0, 12.0, 10.0, 8.0] if args.quick else [20.0, 16.0, 13.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.5, 6.0]
    fits = []
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for (n, m, U) in cases:
        rule = theory.selection_rule(m, n)
        s_pred = rule["rings"][0][2] if rule["rings"] else np.nan
        v = rule["central"]
        print(f"\n(n,m)=({n},{m}), U={U}: predicted central {v}, ring exponent s={s_pred:.3f}, "
              f"slope {-n * s_pred:.3f}")
        rho0 = rhos[0]
        M0 = NonlocalModel(polygon_sample(n, args.h, rho=rho0))
        first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [3.0, 3.5, 2.5, 4.0, 2.0],
                                 verbose=False)
        if first is None:
            print("   no initial solution")
            continue
        s0 = newton.solve_robust(M0, first["w"], first["b"], fixed_power=U)
        if not s0["converged"]:
            print("   fixed-power solve failed at rho0")
            continue

        def build(rho, n=n):
            return NonlocalModel(polygon_sample(n, args.h, rho=rho))

        def solve_fixed_U(model, w0, prev, U=U):
            return newton.solve_robust(model, w0, prev["b"], fixed_power=U)

        scan = [(rho0, M0, s0)] + param_continuation(build, rho0, M0, s0, rhos[1:], solve_fixed_U)
        rows = []
        for rho, M, s in scan:
            row = dg.mode_summary(M, s, nsym=n, m=m)
            row.update(n=n, m=m, rho=rho, eps=theory.index_anisotropy_coefficient(n) * U / rho ** n)
            rows.append(row)
            print(f"   rho={rho:5.2f}: b={s['b']:.4f} W_small={row['W_small']:+.0f} "
                  f"sats={row['n_satellites']} r_s(det)={row['satellite_r_mean']:.4f} "
                  f"r_s(fit)={row['r_s_fit']:.4f}")
        io_utils.write_csv(os.path.join(out, f"scaling_n{n}_m{m}.csv"), rows)
        rr = np.array([r["rho"] for r in rows])
        rs = np.array([r["r_s_fit"] for r in rows])
        ok = np.isfinite(rs) & (rs > 0)
        if ok.sum() >= 3:
            slope = np.polyfit(np.log(rr[ok]), np.log(rs[ok]), 1)[0]
            fits.append(dict(n=n, m=m, U=U, central_pred=v, s_pred=s_pred, slope_pred=-n * s_pred,
                             slope_fit=slope, rho_min=rr[ok].min(), rho_max=rr[ok].max()))
            print(f"   fitted slope {slope:.3f}  vs predicted {-n * s_pred:.3f}")
            ax.loglog(rr[ok], rs[ok], "o-", label=f"n={n}, m={m}: fit {slope:.2f}, theory {-n * s_pred:.2f}")
    ax.set_xlabel("conformal radius rho")
    ax.set_ylabel("satellite radius r_s (harmonic fit)")
    ax.legend(fontsize=7)
    ax.set_title("r_s at fixed power U")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "scaling.png"), dpi=140)
    plt.close(fig)
    io_utils.write_csv(os.path.join(out, "scaling_fits.csv"), fits)
    print("\nE2b done.")


if __name__ == "__main__":
    main()
