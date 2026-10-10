"""
E1b —— 受对称性保护的模式: 基模 (m=0) 与单位涡旋 (m=1) 在正多边形样品中的形变
=============================================================================
理论: m=1 属于 C_n (n>=3) 的复表示 mu=1, 与其共轭 mu=-1 不耦合 -> 中心奇点荷 1 受保护,
但光强与相位被 C_n 调制 (“非正则涡旋”)。一阶微扰: 调制深度 ∝ 各向异性强度 ∝ U (r_ring/rho)^n。

本实验:
  1) 在 poly3..poly6 (保角半径 rho=8) 与圆盘 (R=8) 中求 m=0, m=1 解族 U(b)
     —— 相同 rho 的样品 U(b) 曲线几乎重合 (证明 rho 是正确的“尺寸”度量, 各向同性响应相同)
  2) m=1 涡旋亮环上光强的 cos(n theta) 调制深度 vs U, 不同 n 的标度
  3) 相位非正则度: 绕中心一周 d(theta)/d(phi) 的最大偏离
  4) 每个解: C_n 权重 (mu=1 应 = 1)、中心荷 (应 = 1)

运行:
  python e1b_protected_modes.py --quick   (约 3-5 分钟)
  python e1b_protected_modes.py           (约 30-60 分钟)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import build_sample, finalize_args, robust_family, save_family, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def ring_modulation(model, sol, n, center=(0.0, 0.0)):
    """亮环半径处 |q| 的 cos(n theta) 调制深度 |c_n|/|c_0|, 以及相位梯度的非正则度。"""
    q = dg.field_grid(model, sol["w"])
    rr = dg.ring_radius(q, model.grid, center, rmax=8.0)
    ls, c = dg.angular_spectrum(np.abs(q) + 0j, model.grid, center, rr, lmax=2 * n)
    ls = list(ls)
    mod = abs(c[ls.index(n)]) / abs(c[ls.index(0)])
    th = np.linspace(0, 2 * np.pi, 721)
    vals = dg.interp(q, model.grid, center[0] + 0.6 * rr * np.cos(th), center[1] + 0.6 * rr * np.sin(th))
    ph = np.unwrap(np.angle(vals))
    dphase = np.gradient(ph, th)
    return float(mod), float(np.max(np.abs(dphase - np.mean(dphase)))), float(rr)


def main():
    p = standard_args(__doc__)
    args = finalize_args(p, "E1_symmetry_transfer")
    out, dat = args.res_dir, args.data_dir
    rho = 8.0
    samples = ["poly3", "poly4", "poly5", "poly6", "disk"]
    b_lo, b_hi = (0.5, 4.0) if args.quick else (0.2, 6.0)
    db = 0.5 if args.quick else 0.2
    all_rows = []
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    for kind in samples:
        sample, nsym = build_sample(kind, args.h, rho)
        M = NonlocalModel(sample)
        print("\n" + M.describe())
        for m in (0, 1):
            name = f"e1b_{kind}_m{m}"
            sols, info = robust_family(M, lambda W, m=m: guesses.lg_vortex(M, m, W), b_lo, b_hi, db=db,
                                       b_start=2.0, real=(m == 0), verbose=False)
            print(f"  m={m}: {len(sols)} solutions, b in [{info.get('b_min')}, {info.get('b_max')}]; "
                  f"low end: {info['stop_down']}")
            rows = save_family(M, sols, out, dat, name, nsym=nsym, m=m, plot_bs=[2.0])
            for r, s in zip(rows, sols):
                r.update(sample=kind, m=m, rho=rho)
                if m == 1 and nsym:
                    r["ring_mod"], r["phase_nonuniformity"], _ = ring_modulation(M, s, nsym)
            io_utils.write_csv(os.path.join(out, f"{name}_summary.csv"), rows)
            all_rows += rows
            ls = "-" if m == 0 else "--"
            axs[0].plot([r["b"] for r in rows], [r["U"] for r in rows], ls, label=f"{kind} m={m}")
            if m == 1 and nsym:
                axs[1].loglog([r["U"] for r in rows], [r["ring_mod"] for r in rows], "o-", ms=3,
                              label=f"n={nsym}")
    axs[0].set_xlabel("b")
    axs[0].set_ylabel("U")
    axs[0].legend(fontsize=6)
    axs[0].set_title(f"U(b), all samples rho={rho:g}")
    axs[1].set_xlabel("U")
    axs[1].set_ylabel("cos(n theta) modulation of |q| on ring")
    axs[1].legend(fontsize=7)
    axs[1].set_title("m=1 vortex: C_n modulation")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "e1b_families.png"), dpi=130)
    plt.close(fig)
    io_utils.write_csv(os.path.join(out, "e1b_all.csv"), all_rows)
    print("\nE1b done.")


if __name__ == "__main__":
    main()
