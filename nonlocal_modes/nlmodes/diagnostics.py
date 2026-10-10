"""
diagnostics.py —— 模式诊断量 (power, widths, phase singularities, C_n irreps ...)
================================================================================

本项目最核心的观测量是“相位奇点(涡旋)的拓扑结构”:
  * find_singularities : 在网格小方格(plaquette)上计算相位绕数, 找出所有 ±1 奇点并做亚网格定位
  * winding_on_circle  : 沿任意圆周的相位绕数 (插值, 不需要奇点落在网格上)
  * winding_profile    : 绕数随半径的变化 W(r) —— “嬗变”表现为 W(r) 在 r_s 处跳变
  * angular_spectrum   : 圆周上的角向谐波 c_l(r) —— 判断哪个角动量分量主导
  * cn_weights         : 投影到 C_n 群的不可约表示 (角赝动量 mu), 检验守恒律
  * harmonic_amplitudes: 实函数(如折射率)在圆周上的 cos(k theta)/sin(k theta) 系数 (E1 理论对比)
"""

import json

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates


# ---------------------------------------------------------------------------
# 基本量
# ---------------------------------------------------------------------------
def field_grid(model, w):
    return model.to_grid(np.asarray(w, complex))


def center_of_mass(model, w):
    I = np.abs(w) ** 2
    s = I.sum()
    return float((model.xo * I).sum() / s), float((model.yo * I).sum() / s)


def rms_widths(model, w):
    """积分宽度 W_x = sqrt(∫(x-xc)^2|q|^2 / ∫|q|^2), W_y 同理。"""
    I = np.abs(w) ** 2
    s = I.sum()
    xc, yc = center_of_mass(model, w)
    Wx = np.sqrt(((model.xo - xc) ** 2 * I).sum() / s)
    Wy = np.sqrt(((model.yo - yc) ** 2 * I).sum() / s)
    return float(Wx), float(Wy)


# ---------------------------------------------------------------------------
# 插值 (把网格上的复场插值到任意点)
# ---------------------------------------------------------------------------
def interp(fgrid, grid, xs, ys, order=3):
    """三次样条插值; 网格外取 0。fgrid 可为实数或复数 (Ny,Nx) 数组。"""
    col = (np.asarray(xs) - grid.x[0]) / grid.h
    row = (np.asarray(ys) - grid.y[0]) / grid.h
    coords = np.vstack([np.ravel(row), np.ravel(col)])
    if np.iscomplexobj(fgrid):
        re = map_coordinates(fgrid.real, coords, order=order, mode="constant", cval=0.0)
        im = map_coordinates(fgrid.imag, coords, order=order, mode="constant", cval=0.0)
        out = re + 1j * im
    else:
        out = map_coordinates(fgrid, coords, order=order, mode="constant", cval=0.0)
    return out.reshape(np.shape(xs))


# ---------------------------------------------------------------------------
# 相位奇点
# ---------------------------------------------------------------------------
def _wrap(d):
    return (d + np.pi) % (2 * np.pi) - np.pi


def beam_region(qgrid, grid, rel_threshold=1e-4, smooth=1.5):
    """“光束所在区域”: 高斯平滑后的光强 > rel_threshold*最大值。用于排除远处尾部的数值奇点。"""
    Is = gaussian_filter(np.abs(qgrid) ** 2, sigma=max(smooth / grid.h, 0.5))
    return Is > rel_threshold * Is.max()


def _subcell_zero(f00, f10, f11, f01):
    """在单位方格内求双线性插值复函数的零点 (s,t), 用于奇点亚网格定位。"""
    s, t = 0.5, 0.5
    for _ in range(20):
        f = f00 * (1 - s) * (1 - t) + f10 * s * (1 - t) + f11 * s * t + f01 * (1 - s) * t
        fs = (f10 - f00) * (1 - t) + (f11 - f01) * t
        ft = (f01 - f00) * (1 - s) + (f11 - f10) * s
        # 2x2 实方程组: [Re fs, Re ft; Im fs, Im ft] [ds, dt] = -[Re f, Im f]
        a, bb, c, d = fs.real, ft.real, fs.imag, ft.imag
        det = a * d - bb * c
        if abs(det) < 1e-300:
            break
        ds = (-f.real * d + f.imag * bb) / det
        dt = (-a * f.imag + c * f.real) / det
        s, t = s + ds, t + dt
        if abs(ds) + abs(dt) < 1e-12:
            break
    return float(np.clip(s, 0, 1)), float(np.clip(t, 0, 1))


def _plaquette_charges(q, ambiguous_tol=1e-6):
    """
    每个小方格 a(iy,ix)->b(iy,ix+1)->c(iy+1,ix+1)->d(iy+1,ix) 逆时针的相位绕数。
    返回 (charge 整数数组, ambiguous 布尔数组)。
    若某条边上相位差恰为 ±pi, 绕数没有定义 (ambiguous=True, charge 记 0)。两种来源:
      * 实数场跨过节线 (偶极、多极): 根本没有奇点;
      * 高对称的复数场, 例如电荷 2 的涡旋恰好位于小方格中心 (四角相位依次差 pi)。
    后一种由 find_singularities 用局部加密插值重新判定。
    """
    ph = np.angle(q)
    a, b, c, d = ph[:-1, :-1], ph[:-1, 1:], ph[1:, 1:], ph[1:, :-1]
    steps = [_wrap(b - a), _wrap(c - b), _wrap(d - c), _wrap(a - d)]
    charge = np.rint(sum(steps) / (2 * np.pi)).astype(int)
    ambiguous = np.zeros(charge.shape, bool)
    for st in steps:
        ambiguous |= np.abs(np.abs(st) - np.pi) < ambiguous_tol
    charge[ambiguous] = 0
    return charge, ambiguous


def _locate(q, xs, ys, charge):
    """对 charge 非零的小方格做亚网格定位, 返回 [(x, y, q)]。xs, ys 为该子网格的坐标轴。"""
    out = []
    hx, hy = xs[1] - xs[0], ys[1] - ys[0]
    for iy, ix in zip(*np.nonzero(charge)):
        s_, t_ = _subcell_zero(q[iy, ix], q[iy, ix + 1], q[iy + 1, ix + 1], q[iy + 1, ix])
        out.append((float(xs[ix] + s_ * hx), float(ys[iy] + t_ * hy), int(charge[iy, ix])))
    return out


def find_singularities(qgrid, grid, region=None, exclude_center=None, exclude_radius=0.0,
                       zero_tol=1e-12, refine=8):
    """
    返回单位相位奇点列表 [(x, y, charge), ...]。
    region          : bool 网格, 只在 region 内找 (应排除光学区域以外的点, 那里的 0 是边界条件而非奇点)
    exclude_center/exclude_radius: 不在这个圆内找 (中心区域改用 winding_on_circle 判定, 见下)

    方法: 每个小方格逆时针累加相位差/2pi 得到绕数, 非零者再用双线性插值求零点做亚网格定位。
    重要局限 (混叠): 4 个角点组成的回路无法表示 |荷| >= 2 的奇点 —— 例如中心在方格内的正则 m=3 涡旋
    会被误读为“中心 -1 + 4 个 +1”。所以多重荷只可能出现在对称中心时, 必须排除中心小圆,
    中心荷改用多点圆周绕数 W(r) 判定 (mode_summary 就是这样做的)。
    若某网格点上场严格为 0, 其周围 2x2 方块用三次插值加密后重新检测 (细网格错开非对称的偏移量)。
    实数场 (偶极、多极) 没有相位奇点, 直接返回空列表。
    """
    q = np.asarray(qgrid)
    amp = np.abs(q)
    if np.max(np.abs(q.imag)) <= 1e-12 * amp.max():
        return []
    charge, _ = _plaquette_charges(q)
    reg = np.ones(charge.shape, bool)
    if region is not None:
        reg = region[:-1, :-1] & region[:-1, 1:] & region[1:, 1:] & region[1:, :-1]
    if exclude_center is not None and exclude_radius > 0:
        xc = 0.5 * (grid.x[:-1] + grid.x[1:])
        yc = 0.5 * (grid.y[:-1] + grid.y[1:])
        XC, YC = np.meshgrid(xc, yc)
        reg &= np.hypot(XC - exclude_center[0], YC - exclude_center[1]) >= exclude_radius
    charge[~reg] = 0
    sings = []
    h = grid.h
    Ny, Nx = q.shape
    zeros = amp <= zero_tol * amp.max()
    for iy, ix in zip(*np.nonzero(zeros)):
        if iy < 1 or ix < 1 or iy > Ny - 2 or ix > Nx - 2 or not reg[iy - 1:iy + 1, ix - 1:ix + 1].all():
            continue
        charge[iy - 1:iy + 1, ix - 1:ix + 1] = 0
        sub = 2 * h * (np.arange(2 * refine) + 0.37) / (2 * refine)
        xs, ys = grid.x[ix] - h + sub, grid.y[iy] - h + sub
        XS, YS = np.meshgrid(xs, ys)
        qf = interp(q, grid, XS, YS)
        sings += _locate(qf, xs, ys, _plaquette_charges(qf)[0])
    sings += _locate(q, grid.x, grid.y, charge)
    return sings


def winding_on_circle(qgrid, grid, center, r, npts=None):
    """沿圆周 (center, r) 的相位绕数 (插值得到, 返回浮点数, 正常情况下接近整数)。实数场返回 0。"""
    if np.max(np.abs(np.imag(qgrid))) <= 1e-12 * np.max(np.abs(qgrid)):
        return 0.0
    if npts is None:
        npts = int(max(720, 40 * 2 * np.pi * r / grid.h))
    th = np.linspace(0, 2 * np.pi, npts + 1)
    vals = interp(qgrid, grid, center[0] + r * np.cos(th), center[1] + r * np.sin(th))
    return float(_wrap(np.diff(np.angle(vals))).sum() / (2 * np.pi))


def winding_profile(qgrid, grid, center, radii):
    """W(r): 各半径圆周上的绕数。"""
    return np.array([winding_on_circle(qgrid, grid, center, r) for r in radii])


def angular_spectrum(qgrid, grid, center, r, lmax=12, npts=1024):
    """圆周上 q(r,theta) = sum_l c_l e^{i l theta}, 返回 (l 数组, c_l 数组)。"""
    th = 2 * np.pi * np.arange(npts) / npts
    vals = interp(qgrid, grid, center[0] + r * np.cos(th), center[1] + r * np.sin(th))
    c = np.fft.fft(vals) / npts
    ls = np.arange(-lmax, lmax + 1)
    return ls, c[ls % npts]


def ring_radius(qgrid, grid, center, rmax=None, nr=200):
    """方位平均光强最大处的半径 (涡旋“亮环”半径)。"""
    if rmax is None:
        rmax = 0.5 * min(grid.x[-1] - grid.x[0], grid.y[-1] - grid.y[0])
    rs = np.linspace(0, rmax, nr)
    th = np.linspace(0, 2 * np.pi, 128, endpoint=False)
    # 一次插值所有 (r, theta) 点 (三次样条的预滤波只做一次; 细网格上比逐个半径插值快 ~100 倍)
    v = interp(qgrid, grid, center[0] + np.outer(rs, np.cos(th)), center[1] + np.outer(rs, np.sin(th)))
    Iavg = np.mean(np.abs(v) ** 2, axis=1)
    return float(rs[int(np.argmax(Iavg))])


def estimate_core_radius(qgrid, grid, center, l_out, l_in, radii):
    """
    用中心附近的角向谐波拟合 c_l(r) ≈ a_l r^{|l|}, 估计“嬗变核”半径
        r_s = (|a_{l_in}| / |a_{l_out}|)^{1/(|l_out| - |l_in|)}。
    适用于卫星奇点距中心小于网格间距、无法直接分辨的情形 (亚网格估计)。
    """
    lmax = max(abs(l_out), abs(l_in)) + 1
    a_out, a_in = [], []
    for r in radii:
        ls, c = angular_spectrum(qgrid, grid, center, r, lmax=lmax)
        a_out.append(abs(c[list(ls).index(l_out)]) / r ** abs(l_out))
        a_in.append(abs(c[list(ls).index(l_in)]) / r ** abs(l_in))
    a_out, a_in = np.median(a_out), np.median(a_in)
    if abs(l_out) == abs(l_in) or a_out == 0:
        return np.nan
    return float((a_in / a_out) ** (1.0 / (abs(l_out) - abs(l_in))))


# ---------------------------------------------------------------------------
# C_n 不可约表示投影 (角赝动量)
# ---------------------------------------------------------------------------
def cn_weights(qgrid, grid, n, center=(0.0, 0.0)):
    """
    把场投影到 C_n 的不可约表示: Pi_mu q = (1/n) sum_k e^{i 2 pi mu k/n} R^k q,
    其中 (R^k q)(r,theta) = q(r, theta - 2 pi k/n)。对 q ∝ e^{i m theta} 有 Pi_mu q = delta_{mu, m mod n} q。
    返回 {mu: 权重}, mu 取值于 (-n/2, n/2]。
    """
    Xc, Yc = grid.X - center[0], grid.Y - center[1]
    rots = []
    for k in range(n):
        a = 2 * np.pi * k / n
        xr = center[0] + Xc * np.cos(a) + Yc * np.sin(a)
        yr = center[1] - Xc * np.sin(a) + Yc * np.cos(a)
        rots.append(interp(qgrid, grid, xr, yr) if k else np.asarray(qgrid, complex))
    weights = {}
    tot = 0.0
    for mu in range(n):
        proj = sum(np.exp(2j * np.pi * mu * k / n) * rots[k] for k in range(n)) / n
        wmu = float(np.sum(np.abs(proj) ** 2))
        mu_c = mu if mu <= n // 2 else mu - n
        weights[mu_c] = wmu
        tot += wmu
    return {k: v / tot for k, v in sorted(weights.items())}


def harmonic_amplitudes(fgrid, grid, center, r, kmax=12, npts=1024):
    """实函数 f(r,theta) = a_0 + sum_k [a_k cos k theta + b_k sin k theta], 返回 (a, b) 数组。"""
    th = 2 * np.pi * np.arange(npts) / npts
    vals = interp(np.nan_to_num(fgrid), grid, center[0] + r * np.cos(th), center[1] + r * np.sin(th))
    c = np.fft.rfft(vals) / npts
    a = 2 * c.real
    a[0] /= 2
    b = -2 * c.imag
    return a[:kmax + 1], b[:kmax + 1]


# ---------------------------------------------------------------------------
# 汇总 (写入 CSV 的一行)
# ---------------------------------------------------------------------------
def singularity_structure(sings, center, r_center):
    """把奇点分成“中心”(距中心 < r_center) 与“卫星”, 返回 dict。"""
    central = 0
    sats = []
    for (x, y, q) in sings:
        r = np.hypot(x - center[0], y - center[1])
        if r < r_center:
            central += q
        else:
            sats.append((round(r, 4), round(float(np.arctan2(y - center[1], x - center[0])), 4), q))
    sats.sort()
    return dict(central_charge=central, n_satellites=len(sats),
                satellite_charge=int(sum(s[2] for s in sats)),
                satellite_r_mean=float(np.mean([s[0] for s in sats])) if sats else 0.0,
                satellite_r_min=float(min(s[0] for s in sats)) if sats else 0.0,
                satellite_r_max=float(max(s[0] for s in sats)) if sats else 0.0,
                satellites=json.dumps(sats))


def mode_summary(model, sol, nsym=None, center=(0.0, 0.0), m=None, r_center=None):
    # center=None: 以光束质心为中心 (用于不在样品中心的态, 例如表面态)
    """
    一个定态解的全部标量诊断量 (CSV 的一行)。
    nsym   : 样品的旋转对称阶数 n (若给出则计算 C_n 权重, 以及理论预言的中心荷)
    center : 对称中心 (多边形/圆样品为原点)
    m      : 该解族的拓扑荷标签 (给出时计算嬗变核半径的亚网格拟合值 r_s_fit)

    关键列:
      residual       定态方程残差 (最大模), 检查解的精度
      total_charge   亮环半径处的绕数 (总拓扑荷)
      W_small        半径 0.5h 小圆上的绕数 (插值; 卫星离中心 > 0.5h 时等于真实中心荷)
      central_charge 中心荷 = round(W_small) (卫星环 >= 2h 时可信)
      W_core / unresolved_charge  半径 r_center (默认 1.25h) 圆周绕数; 与 W_small 之差 = 未分辨卫星的荷
      n_satellites / satellites  r_center < r < 1.5 r_ring 处网格检测到的单位奇点 (r, theta, q) 列表
      n_outer_singularities      更外面 (光强指数小的尾部) 的奇点个数, 物理上不重要
      r_s_fit        由中心附近角谐波拟合的嬗变核半径 (可小于网格间距)
      mu_dominant / mu_weight    C_n 不可约表示中占主导者及其权重
    """
    g = model.grid
    w = sol["w"]
    q = field_grid(model, w)
    n_t = sol.get("n_t")
    if n_t is None:
        n_t = model.thermal_index(np.abs(w) ** 2)
    xc, yc = center_of_mass(model, w)
    Wx, Wy = rms_widths(model, w)
    if center is None:
        center = (xc, yc)
    # residual: 定态方程残差的最大模 max|1/2 D w + (V-b) w| (与牛顿收敛判据相同的量);
    # 正常应 < 牛顿判据 (main 1e-9, hq 1e-12); 舍入误差下限约 1e-13~1e-12 (见 docs/numerics.md)
    R = 0.5 * (model.D @ w) + (model.potential(n_t) - sol["b"]) * w
    row = dict(b=sol["b"], U=sol["power"], h=g.h, peak=float(np.abs(w).max()),
               n_max=float(np.max(n_t)), xc=xc, yc=yc, Wx=Wx, Wy=Wy,
               converged=bool(sol.get("converged", True)), iterations=sol.get("iterations", -1),
               residual=float(np.max(np.abs(R))))
    if r_center is None:
        r_center = 1.25 * g.h
    rr = ring_radius(q, g, center, rmax=3 * max(Wx, Wy))
    row["ring_radius"] = rr
    region = beam_region(q, g) & model.opt_mask
    sings = find_singularities(q, g, region, exclude_center=center, exclude_radius=r_center)
    # 卫星只统计亮环附近以内 (r < 1.5 r_ring); 更外面的是指数小尾巴里的零点 (强加热时常见), 单独计数
    r_lim = 1.5 * max(rr, 2 * g.h)
    inner = [sg for sg in sings if np.hypot(sg[0] - center[0], sg[1] - center[1]) < r_lim]
    row.update(singularity_structure(inner, center, r_center))
    row["n_outer_singularities"] = len(sings) - len(inner)
    row["total_charge"] = round(winding_on_circle(q, g, center, max(rr, 2 * g.h)), 3)
    # 中心荷: 多点圆周绕数 (不受方格混叠影响)。注意亚网格尺度不可信: 三次插值在 r < h 处对 r^|m|
    # 型高阶零点 (|m|>=3) 误差可达 O(1)。因此:
    #   W_small (r=0.5h) -> central_charge, 只有当卫星环半径 >= 2h (见 r_s_fit) 时才算“已分辨”;
    #   W_core  (r=r_center) - W_small = 落在 [0.5h, r_center] 内、网格无法分辨的卫星总荷;
    #   W_profile: r = 0.5h, h, 1.5h, 2h, 3h 处的绕数, 供人工判读。
    radii = [0.5 * g.h, g.h, 1.5 * g.h, 2 * g.h, 3 * g.h]
    Wp = [round(winding_on_circle(q, g, center, r), 2) for r in radii]
    row["W_profile"] = json.dumps(Wp)
    row["W_small"] = Wp[0]
    row["W_core"] = round(winding_on_circle(q, g, center, r_center), 3)
    row["central_charge"] = int(round(row["W_small"]))
    row["unresolved_charge"] = int(round(row["W_core"] - row["W_small"]))
    if nsym:
        wts = cn_weights(q, g, nsym, center)
        mu = max(wts, key=wts.get)
        row["mu_dominant"] = mu
        row["mu_weight"] = wts[mu]
    if m is not None:
        from .theory import selection_rule
        v = selection_rule(m, nsym)["central"] if nsym else m
        row["m"] = m
        row["predicted_central"] = v
        row["r_s_fit"] = (estimate_core_radius(q, g, center, m, v, [g.h, 1.5 * g.h, 2 * g.h])
                          if v != m else 0.0)
    return row
