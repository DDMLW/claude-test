"""
E0d —— 网格对称性造成的“假嬗变”有多大? (论文可信度的关键检验)
==============================================================
笛卡尔网格本身具有 C_4 对称性。即使使用各向同性九点格式, 仍残留 O(h^4) 的 C_4 型离散误差。
它会像一个微弱的“正方形样品”那样耦合角谐波 l 与 l±4:
    * m = 3:  m ≡ -1 (mod 4) -> 可能出现网格诱导的“中心 -1 + 4 卫星”, 半径 ∝ eps_grid^(1/2) ∝ h^2
    * m = 4:  m ≡ 0  (mod 4) -> 可能出现网格诱导的“中心 0 + 4 卫星”,  半径 ∝ eps_grid^(1/4) ∝ h
本实验给出“数值假象的底线”, 物理结论必须远高于此底线。

  (a) 圆盘 (物理上 C_inf, 不应有任何嬗变) 中 m = 2,3,4,5 涡旋: 卫星半径 r_s 随 h 的变化
  (b) 旋转检验: 把正多边形相对网格转一个角度 theta0, 物理结构应随样品整体刚性转动、
      卫星半径不变; 网格假象则不随之转动。对 (n,m) = (3,2), (5,3), (6,4) 做 theta0 = 0, 7, 15 度。

运行:
  python e0d_grid_artifacts.py --quick   (约 5 分钟)
  python e0d_grid_artifacts.py           (约 30-60 分钟)
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, pick, pick_closest, robust_family, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils  # noqa: E402
from nlmodes.geometry import disk_sample, polygon_sample  # noqa: E402
from nlmodes.model import NonlocalModel  # noqa: E402


def solve_at(M, m, b_target):
    sols, _ = robust_family(M, lambda W: guesses.lg_vortex(M, m, W), b_target, b_target + 2.0, db=0.25,
                            b_start=b_target + 1.0, verbose=False)
    s = pick_closest(sols, "b", b_target)
    return s if s is not None and abs(s["b"] - b_target) < 1e-6 else None


def main():
    p = standard_args(__doc__)
    p.add_argument("--b", type=float, default=1.0, help="比较所用的传播常数")
    args = finalize_args(p, "E0_validation")
    out = args.res_dir
    rows = []
    hs = pick(args, quick=[0.4, 0.3], main=[0.3, 0.2, 0.15], hq=[0.3, 0.2, 0.15, 0.1])
    print("\n(a) disk control: grid-induced satellites")
    for h in hs:
        M = NonlocalModel(disk_sample(8.0, h))
        for m in ([2, 3, 4] if args.quick else [2, 3, 4, 5]):
            s = solve_at(M, m, args.b)
            if s is None:
                print(f"   h={h} m={m}: no solution")
                continue
            r = dg.mode_summary(M, s, nsym=4, m=m)
            rows.append(dict(part="disk_control", h=h, m=m, b=s["b"], central=r["central_charge"],
                             n_sat=r["n_satellites"], r_s=r["satellite_r_mean"], W_small=r["W_small"]))
            print(f"   h={h} m={m}: central {r['central_charge']:+d}, grid-induced satellites {r['n_satellites']} "
                  f"at r={r['satellite_r_mean']:.4f}")
    print("\n(b) rotation test: polygon rotated against the grid")
    h = pick(args, quick=0.3, main=0.2, hq=0.1)
    for (n, m) in ([(3, 2), (5, 3)] if args.quick else [(3, 2), (5, 3), (6, 4)]):
        for deg in (0.0, 7.0, 15.0):
            M = NonlocalModel(polygon_sample(n, h, rho=8.0, theta0=np.radians(deg)))
            s = solve_at(M, m, args.b)
            if s is None:
                print(f"   n={n} m={m} theta0={deg}: no solution")
                continue
            r = dg.mode_summary(M, s, nsym=None, m=m)
            sats = json.loads(r["satellites"])
            rel = sorted(round(np.degrees(a) - deg, 1) % (360 / n) for (_, a, _) in sats)
            rows.append(dict(part="rotation", n=n, m=m, theta0_deg=deg, h=h, b=s["b"], central=r["central_charge"],
                             n_sat=r["n_satellites"], r_s=r["satellite_r_mean"],
                             sat_angles_rel_mod=json.dumps(rel)))
            print(f"   n={n} m={m} theta0={deg:4.1f}: central {r['central_charge']:+d}, sats {r['n_satellites']}, "
                  f"r_s={r['satellite_r_mean']:.4f}, satellite angles relative to sample (mod 360/n): {rel}")
    io_utils.write_csv(os.path.join(out, "e0d_grid_artifacts.csv"), rows)
    print("\nE0d done.")


if __name__ == "__main__":
    main()
