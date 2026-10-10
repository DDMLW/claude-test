"""
theory.py —— 理论预言 (用于和数值结果自动对比)
===============================================

详细推导见 docs/theory.md。这里只实现可以直接算出数字的公式:

1) 对称性传递 (symmetry transfer):
   正 n 边形样品 (保角半径 rho, 一个顶点在 theta=0 方向), 中心处总能流 U 的 C_n 对称光束
   产生的折射率中, 最低阶的各向异性项为
        dn_aniso(r, theta) ≈ U / (pi n (n+1)) * (r/rho)^n * cos(n theta)
   (index_anisotropy_coefficient / predicted_cos_harmonic)

2) 角赝动量选择定则 + “牛顿多边形”(Newton polygon) 规则:
   电荷 m 的涡旋放进 C_n 样品后, 场只含角谐波 l = m - j n (j = 0,1,2,...), 第 j 个谐波的振幅
   ∝ eps^j (eps 为各向异性强度)。在小 r 处 r^{|l|} eps^j 最大者主导, 由点集 (|l|, j) 的下凸包
   决定: 中心荷 v、卫星奇点环的数目与半径标度 r_ring ∝ eps^{s}。
   (selection_rule)
"""

import numpy as np

from .geometry import polygon_circumradius_from_rho, polygon_integral_In, polygon_rho_from_circumradius  # noqa: F401


def index_anisotropy_coefficient(n):
    """c_n = 1/(pi n (n+1)), 使 dn_aniso = c_n U (r/rho)^n cos(n theta)。"""
    return 1.0 / (np.pi * n * (n + 1))


def predicted_cos_harmonic(n, U, rho, r):
    """理论预言: 折射率在半径 r 圆周上 cos(n theta) 分量的振幅。"""
    return index_anisotropy_coefficient(n) * U * (np.asarray(r) / rho) ** n


def irrep_label(m, n):
    """角赝动量 mu = m mod n, 取值于 (-n/2, n/2]。"""
    mu = m % n
    if mu > n / 2:
        mu -= n
    return int(mu)


def newton_polygon(m, n, jmax=None):
    """
    返回下凸包顶点列表 [(|l|, j, l), ...], 从外层 (l=m, j=0) 到中心主导谐波。
    只考虑 |l| <= |m| 的谐波 l = m - j n (若 m<0 则 l = m + j n)。
    相同 |l| 取 j 最小者 (低阶先主导)。
    """
    if m == 0:
        return [(0, 0, 0)]
    sgn = 1 if m > 0 else -1
    if jmax is None:
        jmax = 2 * abs(m) // n + 2
    pts = {}
    for j in range(jmax + 1):
        l = m - sgn * j * n
        if abs(l) <= abs(m):
            pts.setdefault(abs(l), (abs(l), j, l))
    pts = sorted(pts.values(), key=lambda p: -p[0])          # |l| 从大到小
    hull = [pts[0]]
    for p in pts[1:]:
        # 维护下凸包: 斜率 (j 增量)/(|l| 减量) 必须递增
        while len(hull) >= 2:
            a, bpt = hull[-2], hull[-1]
            s1 = (bpt[1] - a[1]) / (a[0] - bpt[0])
            s2 = (p[1] - bpt[1]) / (bpt[0] - p[0])
            if s2 <= s1 + 1e-12:
                hull.pop()
            else:
                break
        hull.append(p)
    # 只保留沿 |l| 减小时 j 的“有效”顶点 (凸包下沿到 |l| 最小点为止)
    return hull


def selection_rule(m, n):
    """
    电荷 m 的涡旋在 C_n 样品中的结构预言 (弱各向异性极限)。
    返回 dict:
        mu            角赝动量 (守恒量)
        central       中心奇点的荷 v
        rings         [(卫星个数, 每个卫星的荷, 半径标度指数 s), ...]  r_ring ∝ eps^s, 由外到内
        degenerate    n | 2m 时 +m 与 -m 被各向异性直接耦合 (C_nv 下涡旋需非线性锁相, 预期存在阈值)
        n_satellites  卫星总数
    n = 0 或 None 表示圆对称 (C_inf): 没有嬗变。
    """
    if not n:
        return dict(mu=m, central=m, rings=[], degenerate=False, n_satellites=0)
    hull = newton_polygon(m, n)
    rings = []
    for a, bpt in zip(hull[:-1], hull[1:]):
        s = (bpt[1] - a[1]) / (a[0] - bpt[0])
        dq = a[2] - bpt[2]
        rings.append((abs(dq), int(np.sign(dq)), float(s)))
    degenerate = (m != 0) and ((2 * m) % n == 0)
    return dict(mu=irrep_label(m, n), central=int(hull[-1][2]), rings=rings,
                degenerate=bool(degenerate), n_satellites=int(sum(r[0] for r in rings)))


def selection_table(ms=range(1, 7), ns=(2, 3, 4, 5, 6)):
    """打印/返回整张“周期表”预言 (用于 E2 结果对照)。"""
    rows = []
    for n in ns:
        for m in ms:
            r = selection_rule(m, n)
            rows.append(dict(n=n, m=m, mu=r["mu"], central=r["central"],
                             n_satellites=r["n_satellites"],
                             rings=";".join(f"{k}x({q:+d}) ~eps^{s:.3g}" for k, q, s in r["rings"]),
                             degenerate=r["degenerate"]))
    return rows
