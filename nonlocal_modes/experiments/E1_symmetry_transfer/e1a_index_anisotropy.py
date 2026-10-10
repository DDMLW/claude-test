"""
E1a —— “对称性传递”定律: 均匀热介质的折射率继承样品边界的 C_n 对称性
=====================================================================
理论 (docs/theory.md 第 2 节): 正 n 边形样品 (保角半径 rho, 一个顶点指向 theta=0) 中心处
任意 C_n 对称、总能流为 U 的光束所产生的折射率, 其最低阶各向异性项为

        dn_n(r) cos(n theta),   dn_n(r) = U / (pi n (n+1)) * (r / rho)^n

关键性质: 系数只依赖 U 与 rho, 与光束形状无关 (普适性); 次级谐波只有 2n, 3n, ...

本实验数值验证:
  1) 幂律与系数: 对 n = 3,4,5,6, 不同 rho, 拟合 dn_n(r) 与理论比较 (log-log 斜率 = n)
  2) 普适性: 同一 U 下改变光束形状 (高斯宽度、环形光束), 系数不变
  3) 线性: dn_n 正比于 U
  4) 选择性: 只有 k = 0, n, 2n, ... 的角谐波非零 (其余为数值噪声)
  5) 矩形 (C_2) 样品: 没有闭式解, 给出数值系数 c_2(纵横比) 供后续使用

运行:
  python e1a_index_anisotropy.py --quick   (约 1 分钟)
  python e1a_index_anisotropy.py           (约 5-10 分钟)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import io_utils, theory  # noqa: E402
from nlmodes.geometry import polygon_sample, rectangle_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def gaussian(M, w, U):
    I = np.exp(-(M.xo ** 2 + M.yo ** 2) / w ** 2)
    return I * U / (M.h ** 2 * I.sum())


def ring(M, R0, w, U):
    r = np.hypot(M.xo, M.yo)
    I = np.exp(-((r - R0) / w) ** 2)
    return I * U / (M.h ** 2 * I.sum())


def harmonics(M, I, r, kmax):
    ng = M.index_grid(M.thermal_index(I))
    a, b = dg.harmonic_amplitudes(ng, M.grid, (0, 0), r, kmax=kmax)
    return a, b


def main():
    p = standard_args(__doc__, default_h=0.2, quick_h=0.3)
    args = finalize_args(p, "E1_symmetry_transfer")
    out = args.res_dir
    h = args.h
    rhos = [8.0, 12.0] if args.quick else [6.0, 8.0, 10.0, 12.0, 15.0]
    radii = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0])
    U = 10.0
    rows = []

    # 1) power law & coefficient
    print("\n[1] dn_n(r) vs theory")
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, n in zip(axs, (3, 4, 5, 6)):
        for rho in rhos:
            M = NonlocalModel(polygon_sample(n, h, rho=rho))
            I = gaussian(M, 1.5, U)
            rr = radii[radii < 0.6 * rho]
            vals = []
            for r in rr:
                a, b = harmonics(M, I, r, 3 * n)
                pred = theory.predicted_cos_harmonic(n, U, rho, r)
                others = [abs(a[k]) for k in range(1, 3 * n + 1) if k % n] + [np.max(np.abs(b))]
                rows.append(dict(part="power_law", n=n, rho=rho, U=U, beam="gauss w=1.5", r=r,
                                 dn_n=a[n], theory=pred, ratio=a[n] / pred, dn_2n=a[2 * n],
                                 max_forbidden=max(others)))
                vals.append(a[n])
            ax.loglog(rr, vals, "o", label=f"rho={rho:g}")
            ax.loglog(rr, theory.predicted_cos_harmonic(n, U, rho, rr), "k-", lw=0.7)
            print(f"   n={n} rho={rho:5.1f}: ratio(num/theory) at r={rr[0]:g}..{rr[-1]:g}: "
                  + " ".join(f"{v / theory.predicted_cos_harmonic(n, U, rho, r):.4f}" for v, r in zip(vals, rr)))
        ax.set_title(f"n={n}: dn_n(r), lines = U/(pi n(n+1)) (r/rho)^n")
        ax.set_xlabel("r")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "e1a_power_law.png"), dpi=130)
    plt.close(fig)

    # 2) universality (beam shape) & 3) linearity in U
    print("\n[2-3] universality and linearity (n=3, rho=10, r=3)")
    M = NonlocalModel(polygon_sample(3, h, rho=10.0))
    pred = theory.predicted_cos_harmonic(3, U, 10.0, 3.0)
    for name, I in [("gauss w=1", gaussian(M, 1.0, U)), ("gauss w=2", gaussian(M, 2.0, U)),
                    ("gauss w=3", gaussian(M, 3.0, U)), ("ring R0=2 w=0.7", ring(M, 2.0, 0.7, U)),
                    ("ring R0=4 w=1", ring(M, 4.0, 1.0, U))]:
        a, _ = harmonics(M, I, 3.0, 6)
        rows.append(dict(part="universality", n=3, rho=10.0, U=U, beam=name, r=3.0, dn_n=a[3], theory=pred,
                         ratio=a[3] / pred))
        print(f"   {name:16s}: ratio = {a[3] / pred:.4f}")
    for Ux in (2.0, 5.0, 20.0, 40.0):
        a, _ = harmonics(M, gaussian(M, 1.5, Ux), 3.0, 6)
        p_ = theory.predicted_cos_harmonic(3, Ux, 10.0, 3.0)
        rows.append(dict(part="linearity", n=3, rho=10.0, U=Ux, beam="gauss w=1.5", r=3.0, dn_n=a[3], theory=p_,
                         ratio=a[3] / p_))
        print(f"   U={Ux:5.1f}: ratio = {a[3] / p_:.4f}")

    # 5) rectangle C_2 coefficient: dn_2(r) = c_2(aspect) * U * (r/rho_eff)^2 ; report c2 = dn_2 / (U r^2)
    print("\n[5] rectangles (C_2): numerical coefficient dn_2/(U r^2) vs aspect ratio (long side 40)")
    for d in ([20, 40] if args.quick else [10, 15, 20, 25, 30, 35, 40]):
        Mr = NonlocalModel(rectangle_sample(d, 40, max(h, 0.25)))
        a, _ = harmonics(Mr, gaussian(Mr, 1.5, U), 2.0, 4)
        rows.append(dict(part="rectangle", n=2, d_eta=d, U=U, r=2.0, dn_n=a[2], coef=a[2] / (U * 4.0)))
        print(f"   {d:g}x40: dn_2(r=2) = {a[2]:+.4e}  ->  coefficient = {a[2] / (U * 4.0):+.4e}")

    io_utils.write_csv(os.path.join(out, "e1a_index_anisotropy.csv"), rows)
    print("\nE1a done.")


if __name__ == "__main__":
    main()
