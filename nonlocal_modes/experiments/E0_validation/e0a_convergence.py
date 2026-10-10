"""
E0a —— 数值方法的收敛性与可靠性检验
====================================
回答一个问题: “我们看到的效应是不是数值误差?” 论文的 Methods/Supplement 需要这些数据。

五个小检验:
  1) 热传导求解精度: 圆盘 + 高斯热源, 与解析解比较, 误差随 h 二阶下降
  2) 对称性传递系数: 正多边形中 cos(n theta) 分量与理论 U/(pi n (n+1)) (r/rho)^n 之比 -> 1
  3) 双涡旋 (Ye et al. OL 2010) 的奇点间距 D_s(b=2) 随 h 收敛
  4) 网格各向异性: 圆盘中 m=2 涡旋的 |c_{-2}/c_2| (纯网格效应), 五点 vs 九点格式
  5) 传播步长: 定态解传播 xi=20 后的偏离度随 dz 二阶下降

运行:
  python e0a_convergence.py --quick      (约 2 分钟)
  python e0a_convergence.py              (约 15-30 分钟)
"""

import os
import sys

import numpy as np
from scipy.special import exp1

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, initial_solution, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, propagate, theory  # noqa: E402
from nlmodes.geometry import disk_sample, polygon_sample, rectangle_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def poisson_error(h, stencil, R=8.0, w=1.0, U=5.0):
    M = NonlocalModel(disk_sample(R, h), stencil=stencil)
    I = U / (np.pi * w ** 2) * np.exp(-(M.xo ** 2 + M.yo ** 2) / w ** 2)
    n = M.potential(M.thermal_index(I))
    r = np.hypot(M.xo, M.yo)
    exact = U / (2 * np.pi) * (np.log(R / r) - 0.5 * exp1(r ** 2 / w ** 2) + 0.5 * exp1(R ** 2 / w ** 2))
    return float(np.max(np.abs(n - exact)))


def transfer_ratio(n, h, rho=10.0, U=10.0, r=3.0):
    M = NonlocalModel(polygon_sample(n, h, rho=rho))
    I = np.exp(-(M.xo ** 2 + M.yo ** 2) / 2.0)
    I *= U / M.power(np.sqrt(I))
    a, _ = dg.harmonic_amplitudes(M.index_grid(M.thermal_index(I)), M.grid, (0, 0), r, kmax=n)
    return float(a[n] / theory.predicted_cos_harmonic(n, U, rho, r))


def twin_vortex_Ds(h, b_target=2.0):
    M = NonlocalModel(rectangle_sample(20, 40, h))
    first = initial_solution(M, lambda W: guesses.lg_vortex(M, 2, W), [4.0, 4.5, 3.5], verbose=False)
    sols, _ = newton.continuation(M, first["w"], first["b"], b_target, -0.25, verbose=False)
    row = dg.mode_summary(M, sols[-1], nsym=2, m=2)
    rs = row["satellite_r_mean"] if row["n_satellites"] == 2 else row["r_s_fit"]
    return 2 * rs, sols[-1]["power"]


def grid_artifact(h, stencil, b=3.0):
    M = NonlocalModel(disk_sample(10, h), stencil=stencil)
    sol = initial_solution(M, lambda W: guesses.lg_vortex(M, 2, W), [b, b - 0.5, b + 0.5], verbose=False)
    q = dg.field_grid(M, sol["w"])
    rr = dg.ring_radius(q, M.grid, (0, 0), rmax=8)
    ls, c = dg.angular_spectrum(q, M.grid, (0, 0), rr, lmax=6)
    ls = list(ls)
    return float(abs(c[ls.index(-2)]) / abs(c[ls.index(2)]))


def main():
    p = standard_args(__doc__)
    args = finalize_args(p, "E0_validation")
    out = args.res_dir
    rows = []

    # 1) Poisson
    hs = [0.4, 0.2, 0.1] if args.quick else [0.4, 0.3, 0.2, 0.15, 0.1]
    print("\n[1] Poisson solver error vs h")
    errs = {5: [], 9: []}
    for h in hs:
        for st in (5, 9):
            e = poisson_error(h, st)
            errs[st].append(e)
            rows.append(dict(test="poisson_disk", h=h, stencil=st, value=e))
            print(f"   h={h:5.3f} stencil={st}  max error = {e:.3e}")

    # 2) symmetry transfer coefficient
    print("\n[2] symmetry-transfer ratio (numerics / theory), should -> 1")
    hs2 = [0.4, 0.2] if args.quick else [0.4, 0.3, 0.2, 0.15]
    for n in (3, 4, 5, 6):
        for h in hs2:
            rt = transfer_ratio(n, h)
            rows.append(dict(test=f"transfer_ratio_n{n}", h=h, stencil=9, value=rt))
            print(f"   n={n} h={h:5.3f}  ratio = {rt:.6f}")

    # 3) twin vortex D_s convergence
    print("\n[3] twin-vortex singularity separation D_s(b=2) in 20x40 rectangle vs h")
    hs3 = [0.4, 0.25] if args.quick else [0.3, 0.25, 0.2, 0.15]
    Ds = []
    for h in hs3:
        d, U = twin_vortex_Ds(h)
        Ds.append(d)
        rows.append(dict(test="twin_vortex_Ds_b2", h=h, stencil=9, value=d, U=U))
        print(f"   h={h:5.3f}  D_s = {d:.5f}   U = {U:.5f}")

    # 4) grid anisotropy artefact
    print("\n[4] grid-induced |c_-2/c_2| for m=2 vortex in a disk (pure numerical artefact)")
    hs4 = [0.4, 0.3] if args.quick else [0.3, 0.2, 0.15]
    for h in hs4:
        for st in (5, 9):
            a = grid_artifact(h, st)
            rows.append(dict(test="grid_artifact_c-2_over_c2", h=h, stencil=st, value=a))
            print(f"   h={h:5.3f} stencil={st}  |c_-2/c_2| = {a:.3e}")

    # 5) propagation time step
    print("\n[5] propagation of a stationary m=1 vortex (triangle rho=8): deviation at xi=20 vs dz")
    M = NonlocalModel(polygon_sample(3, 0.4 if args.quick else 0.25, rho=8.0))
    sol = initial_solution(M, lambda W: guesses.lg_vortex(M, 1, W), [2.0, 2.5, 1.5], verbose=False)
    zmax = 10.0 if args.quick else 20.0
    for dz in ([0.05, 0.025] if args.quick else [0.1, 0.05, 0.025, 0.0125]):
        res = propagate.propagate(M, sol["w"], zmax, dz=dz, w_ref=sol["w"], verbose=False)
        rows.append(dict(test="propagation_deviation", h=M.h, stencil=9, dz=dz,
                         value=float(res["deviation"][-1]),
                         power_drift=float(np.ptp(res["power"]) / res["power"][0])))
        print(f"   dz={dz:7.4f}  deviation(xi={zmax:g}) = {res['deviation'][-1]:.3e}   "
              f"relative power drift = {np.ptp(res['power']) / res['power'][0]:.1e}")

    io_utils.write_csv(os.path.join(out, "convergence_table.csv"), rows)

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    for st, mk in ((5, "s--"), (9, "o-")):
        axs[0].loglog(hs, errs[st], mk, label=f"{st}-point")
    axs[0].loglog(hs, errs[9][-1] * (np.array(hs) / hs[-1]) ** 2, "k:", label="slope 2")
    axs[0].set_xlabel("h")
    axs[0].set_ylabel("max |n - n_exact|")
    axs[0].legend()
    axs[0].set_title("Poisson solver (disk)")
    axs[1].plot(hs3, Ds, "o-")
    axs[1].set_xlabel("h")
    axs[1].set_ylabel("D_s (b=2)")
    axs[1].set_title("twin vortex, 20x40 rectangle")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "convergence.png"), dpi=130)
    print("\nE0a done.")


if __name__ == "__main__":
    main()
