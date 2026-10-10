"""
E5a —— 可重构的热边界调控: 用边界温度多极图样代替样品形状来“设计”涡旋的拓扑结构
================================================================================
思想: 圆形样品 (半径 R) 的边界温度设为 T_b(theta) = Theta * cos(n (theta - theta0))
(实验上可用沿圆周分布的珀尔帖/加热片实现)。热传导方程线性 -> 样品内多出一个与光功率无关的
调和折射率  h(r,theta) = Theta (r/R)^n cos(n(theta - theta0))  —— 纯 C_n 各向异性, 强度连续可调,
符号可反转, 取向 theta0 可旋转。对比: 多边形样品的各向异性 ∝ U/rho^n, 与功率绑定。

预言:
    * 嬗变规则与多边形样品相同 (同为 C_n 对称), 但卫星环半径 r_s ∝ (Theta/R^n)^s 可由 Theta 连续调节
    * Theta 反号 -> 卫星环整体旋转 pi/n
    * Theta -> 0 恢复正则涡旋 (可逆)
    * 调和函数无内部极大值 (极值原理), 外加图样本身不能导引光 —— 只有和自陷效应结合才产生定态

本实验: 对 (n,m) = (3,2), (4,3), (3,3), (5,3), (6,4), (2,2): 固定能流 U, 沿 Theta 做参数延拓 (正负两个方向),
记录 r_s(Theta) (检测值与拟合值)、第一颗卫星的方位角、中心荷。

输出: heating_n<n>_m<m>.csv, heating_rs.png (r_s vs |Theta| 对数图, 带理论斜率), heating_modes_*.png

运行:
  python e5a_boundary_heating.py --quick    (约 5 分钟)
  python e5a_boundary_heating.py            (约 30-60 分钟)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, initial_solution, param_continuation, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting, theory  # noqa: E402
from nlmodes.geometry import disk_sample, harmonic_multipole  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# (n, m, 固定能流 U, 正式运行的最大边界温度幅度 Theta_max)
# n=5,6 在多边形样品 (rho=8) 中物理各向异性太弱, 会被网格的 C_4 误差掩盖 (见 E0d);
# 用加热图样可以把各向异性调得足够强, 是展示 n=5,6 嬗变规则的主要途径。
# n=2 (四极加热) 的中心是外加势的鞍点, Theta 过大会把光束推离中心, 所以上限取小。
CASES = [(3, 2, 15.0, 4.0), (4, 3, 20.0, 4.0), (3, 3, 20.0, 4.0), (5, 3, 20.0, 4.0), (6, 4, 25.0, 4.0),
         (2, 2, 15.0, 1.0)]


def heated_disk(R, h, n, Theta, theta0=0.0):
    s = disk_sample(R, h)
    if Theta != 0:
        s.set_boundary_temperature(harmonic_multipole(n, Theta, R, theta0))
    return NonlocalModel(s)


def first_satellite_angle(sats_json):
    import json
    sats = json.loads(sats_json)
    if not sats:
        return np.nan
    return float(min(sats, key=lambda s: s[0])[1])


def main():
    p = standard_args(__doc__)
    p.add_argument("--R", type=float, default=10.0, help="圆盘半径")
    p.add_argument("--theta_scale", type=float, default=1.0, help="把所有 Theta_max 乘以此因子")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    R = args.R
    nst = 4 if args.quick else 10
    cases = CASES[:2] if args.quick else CASES
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for (n, m, U, tmax) in cases:
        tmax = (1.0 if args.quick else tmax) * args.theta_scale
        thetas = list(np.round(np.geomspace(tmax / 2 ** (nst - 1), tmax, nst), 6))
        rule = theory.selection_rule(m, n)
        s_pred = rule["rings"][0][2] if rule["rings"] else np.nan
        print(f"\n(n,m)=({n},{m}), U={U}: predicted central {rule['central']}, r_s ~ Theta^{s_pred:.3f}")
        M0 = heated_disk(R, args.h, n, 0.0)
        first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [3.0, 3.5, 2.5, 4.0, 2.0],
                                 verbose=False)
        if first is None:
            print("   no initial solution")
            continue
        s0 = newton.solve_robust(M0, first["w"], first["b"], fixed_power=U)

        def solve_fixed_U(model, w0, prev, U=U):
            return newton.solve_robust(model, w0, prev["b"], fixed_power=U)

        rows = []
        for sign in (+1, -1):
            scan = param_continuation(lambda T, n=n: heated_disk(R, args.h, n, T), 0.0, M0, s0,
                                      [sign * t for t in thetas], solve_fixed_U)
            for T, M, s in scan:
                row = dg.mode_summary(M, s, nsym=n, m=m)
                row.update(n=n, m=m, Theta=T, sat_angle=first_satellite_angle(row["satellites"]))
                rows.append(row)
                print(f"   Theta={T:+8.4f}: b={s['b']:.4f} W_small={row['W_small']:+.0f} sats={row['n_satellites']} "
                      f"r_s(det)={row['satellite_r_mean']:.4f} r_s(fit)={row['r_s_fit']:.4f} "
                      f"angle={row['sat_angle']:+.3f}")
                if abs(T) == tmax:
                    plotting.plot_mode(M, s, os.path.join(out, f"heating_mode_n{n}_m{m}_T{T:+g}.png"),
                                       title=f"disk R={R:g}, T_b = {T:+g} cos({n} theta), m={m}, U={U:g}")
        rows.sort(key=lambda r: r["Theta"])
        io_utils.write_csv(os.path.join(out, f"heating_n{n}_m{m}.csv"), rows)
        pos = [r for r in rows if r["Theta"] > 0 and r["r_s_fit"] > 0]
        if len(pos) >= 3:
            T_ = np.array([r["Theta"] for r in pos])
            rs = np.array([r["r_s_fit"] for r in pos])
            slope = np.polyfit(np.log(T_), np.log(rs), 1)[0]
            ax.loglog(T_, rs, "o-", label=f"n={n},m={m}: slope {slope:.2f} (theory {s_pred:.2f})")
            print(f"   fitted slope d ln r_s / d ln Theta = {slope:.3f}  (theory {s_pred:.3f})")
    ax.set_xlabel("boundary temperature amplitude Theta")
    ax.set_ylabel("satellite radius r_s (harmonic fit)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "heating_rs.png"), dpi=140)
    plt.close(fig)
    print("\nE5a done.")


if __name__ == "__main__":
    main()
