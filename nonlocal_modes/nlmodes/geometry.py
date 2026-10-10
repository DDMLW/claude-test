"""
geometry.py —— 计算网格与样品几何 (grid & sample cross-section)
=================================================================

物理模型（归一化单位，与导师论文 Ye et al. PRA 2008 / OL 2010 相同）：

    i dq/dxi = -1/2 (q_xx + q_yy) - [n_lin(x,y) + n(x,y)] q          (光场, 在光学区域内)
    div( kappa grad n ) = -|q|^2                                       (热致折射率, 在热介质内)

本文件只负责“几何”：
  * Grid   : 关于原点对称的“格心”均匀网格, 节点 x_i = (i + 1/2) h (i = -Mx..Mx-1)。
             原点与坐标轴都不是网格点/网格线: 对称解的中心奇点、以及落在镜面对称轴上的卫星奇点
             都位于网格小方格内部, 相位绕数检测不会因为“奇点恰好压在网格点上”而失效。
             网格仍保持关于 x 轴、y 轴、对角线的镜面对称和 C_4 旋转对称。
  * Sample : 一块样品 = 网格 + 每个网格点的“热学类型” + 光学掩膜 + 线性折射率背景。

热学类型 (th_type) 约定：
    MEDIUM    = 0  热非线性介质(温度/折射率未知量, 吸收光并产生热)
    SINK      = 1  热沉/恒温边界(Dirichlet 边界, 其值由 th_bc 给出, 默认 0)
    INSULATOR = 2  绝热材料(零热流, Neumann 边界; 例如导师 PRA 2008 中的线性介质一侧)

光学掩膜 (opt_mask)：True 的点上光场是未知量；掩膜外光场恒为 0 (Dirichlet)。

边界处理：
  * 光场: 未知点 = 落在形状内部的网格点 (阶梯近似)。光束远离边界, 这一近似对光场影响可忽略。
  * 热场: 用水平集 (levelset) 给出真实边界, 采用 Gibou 等人的对称二阶格式 (见 operators.py),
    使热响应(格林函数)的边界误差为 O(h^2), 不会人为破坏样品的 C_n 对称性。
  论文级结果仍需做网格收敛性检验 (见 experiments/E0_validation)。
"""

import numpy as np

MEDIUM, SINK, INSULATOR = 0, 1, 2


class Grid:
    """
    格心均匀网格, 覆盖 [-half_x, half_x] x [-half_y, half_y] (向上取整到 h 的整数倍)。
    数组形状统一为 (Ny, Nx)，第 0 维是 y，第 1 维是 x (与 matplotlib imshow 一致)。
    """

    def __init__(self, half_x, half_y, h):
        self.h = float(h)
        self.Mx = int(np.ceil(half_x / h - 1e-9))
        self.My = int(np.ceil(half_y / h - 1e-9))
        self.x = self.h * (np.arange(-self.Mx, self.Mx) + 0.5)
        self.y = self.h * (np.arange(-self.My, self.My) + 0.5)
        self.X, self.Y = np.meshgrid(self.x, self.y)          # (Ny, Nx)
        self.shape = self.X.shape
        self.R = np.hypot(self.X, self.Y)
        self.TH = np.arctan2(self.Y, self.X)

    @property
    def extent(self):
        """给 imshow 用的 extent (单元中心对齐)。"""
        h = self.h
        return [self.x[0] - h / 2, self.x[-1] + h / 2, self.y[0] - h / 2, self.y[-1] + h / 2]


class Sample:
    """
    一块热非局域样品的全部几何/材料信息。

    属性
    ----
    grid     : Grid
    th_type  : int 数组 (Ny,Nx), 取值 MEDIUM / SINK / INSULATOR
    th_bc    : float 数组 (Ny,Nx), SINK 点上的边界温度(以归一化折射率为单位), 默认 0
    kappa    : float 数组 (Ny,Nx), 相对热导率, 默认 1
    opt_mask : bool 数组 (Ny,Nx), 光场未知量所在的点
    n_lin    : float 数组 (Ny,Nx), 线性折射率背景(例如界面问题中的 n_d = -100)
    name     : 字符串
    params   : dict, 记录构造参数(写入结果文件, 便于复现)
    """

    def __init__(self, grid, th_type, opt_mask, n_lin=None, th_bc=None, kappa=None,
                 name="sample", params=None, levelset=None, bc_func=None):
        self.grid = grid
        # levelset(X, Y): 热介质区域的水平集函数 (内部 < 0, 外部 > 0)。给出时, 热传导算子用
        #   Gibou et al. (J. Comput. Phys. 176, 205 (2002)) 的对称二阶边界格式: 沿网格线线性插值
        #   水平集得到真实边界位置, 只修改对角元 -> 矩阵仍对称, 边界误差 O(h^2) (阶梯近似为 O(h))。
        # bc_func(x, y): 真实边界点上的温度 (给出时优先于 th_bc)。
        self.levelset = levelset
        self.bc_func = bc_func
        self.th_type = np.asarray(th_type, dtype=np.int8)
        self.opt_mask = np.asarray(opt_mask, dtype=bool)
        self.n_lin = np.zeros(grid.shape) if n_lin is None else np.asarray(n_lin, float)
        self.th_bc = np.zeros(grid.shape) if th_bc is None else np.asarray(th_bc, float)
        self.kappa = np.ones(grid.shape) if kappa is None else np.asarray(kappa, float)
        self.name = name
        self.params = dict(params or {})
        self._check()

    def _check(self):
        med = self.th_type == MEDIUM
        if not med.any():
            raise ValueError("样品中没有热介质点 (no MEDIUM nodes)")
        # 热介质点不能贴着网格最外圈, 否则边界条件无定义
        edge = np.zeros(self.grid.shape, bool)
        edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
        if (med & edge).any():
            raise ValueError("热介质点触到了网格边缘, 请增大网格 margin")
        if (self.opt_mask & edge).any():
            raise ValueError("光学未知点触到了网格边缘, 请增大网格 margin")

    # ---- 便捷属性 ----
    @property
    def medium(self):
        return self.th_type == MEDIUM

    def set_boundary_temperature(self, func):
        """
        设置热沉(SINK)点上的温度 = func(X, Y)。
        推荐 func 取一个“调和函数”(Laplace 方程的解, 例如 Theta*(r/R)^n*cos(n*theta)),
        这样阶梯边界上的取值与真实圆周边界上的取值一致, 内部解也就是该调和函数本身。
        """
        sink = self.th_type == SINK
        vals = func(self.grid.X, self.grid.Y)
        self.th_bc = np.where(sink, vals, 0.0)
        self.bc_func = func
        self.params["boundary_temperature"] = getattr(func, "description", "custom")
        return self


# ---------------------------------------------------------------------------
# 形状判定函数 (返回 bool 数组: 是否严格位于形状内部)
# ---------------------------------------------------------------------------
_EPS = 1e-9


def inside_rectangle(X, Y, Lx, Ly):
    return (np.abs(X) < Lx / 2 - _EPS) & (np.abs(Y) < Ly / 2 - _EPS)


def inside_regular_polygon(X, Y, n, R_circ, theta0=0.0):
    """
    正 n 边形, 外接圆半径 R_circ, 一个顶点位于极角 theta0 方向。
    各边外法线方向为 theta0 + pi/n + 2*pi*k/n, 内切圆半径(边心距) = R_circ*cos(pi/n)。
    """
    r_in = R_circ * np.cos(np.pi / n)
    proj_max = np.full(X.shape, -np.inf)
    for k in range(n):
        phi = theta0 + np.pi / n + 2 * np.pi * k / n
        proj_max = np.maximum(proj_max, X * np.cos(phi) + Y * np.sin(phi))
    return proj_max < r_in - _EPS


def inside_disk(X, Y, R):
    return np.hypot(X, Y) < R - _EPS


def inside_ellipse(X, Y, a, b):
    return (X / a) ** 2 + (Y / b) ** 2 < 1.0 - _EPS


# 水平集函数 (内部 < 0, 外部 > 0; 零等值线就是真实边界)
def levelset_rectangle(Lx, Ly):
    return lambda X, Y: np.maximum(np.abs(X) - Lx / 2, np.abs(Y) - Ly / 2)


def levelset_regular_polygon(n, R_circ, theta0=0.0):
    r_in = R_circ * np.cos(np.pi / n)

    def phi(X, Y):
        out = np.full(np.shape(X), -np.inf)
        for k in range(n):
            a = theta0 + np.pi / n + 2 * np.pi * k / n
            out = np.maximum(out, X * np.cos(a) + Y * np.sin(a))
        return out - r_in
    return phi


def levelset_disk(R):
    return lambda X, Y: np.hypot(X, Y) - R


def levelset_ellipse(a, b):
    return lambda X, Y: (np.sqrt((X / a) ** 2 + (Y / b) ** 2) - 1.0) * min(a, b)


# ---------------------------------------------------------------------------
# 正多边形的保角半径 (conformal radius), 理论部分会用到, 见 docs/theory.md
# ---------------------------------------------------------------------------
def polygon_integral_In(n):
    """I_n = int_0^1 (1-s^n)^(-2/n) ds = Gamma(1/n) Gamma(1-2/n) / (n Gamma(1-1/n)),  n>=3。"""
    from math import gamma
    if n < 3:
        raise ValueError("正多边形要求 n >= 3")
    return gamma(1.0 / n) * gamma(1.0 - 2.0 / n) / (n * gamma(1.0 - 1.0 / n))


def polygon_rho_from_circumradius(n, R_circ):
    """正 n 边形(外接圆半径 R_circ)关于中心的保角半径 rho = R_circ / I_n。"""
    return R_circ / polygon_integral_In(n)


def polygon_circumradius_from_rho(n, rho):
    return rho * polygon_integral_In(n)


# ---------------------------------------------------------------------------
# 样品工厂函数
# ---------------------------------------------------------------------------
def _closed_shape_sample(grid, inside, name, params, levelset, n_lin=None):
    """内部点 = 热介质 = 光学未知点; 其余点 = 热沉(温度 0), 光场为 0。"""
    th_type = np.where(inside, MEDIUM, SINK).astype(np.int8)
    return Sample(grid, th_type, inside.copy(), n_lin=n_lin, name=name, params=params,
                  levelset=levelset)


def rectangle_sample(Lx, Ly, h, margin=3):
    """
    矩形样品 [-Lx/2, Lx/2] x [-Ly/2, Ly/2], 四边恒温 (n = q = 0), 对应 Ye et al. OL 2010。
    热边界的真实位置由水平集给出 (Gibou 格式), 不要求 L/2 是 h 的整数倍。
    """
    grid = Grid(Lx / 2 + margin * h, Ly / 2 + margin * h, h)
    inside = inside_rectangle(grid.X, grid.Y, Lx, Ly)
    return _closed_shape_sample(grid, inside, f"rect_{Lx:g}x{Ly:g}",
                                dict(shape="rectangle", Lx=Lx, Ly=Ly, h=h),
                                levelset_rectangle(Lx, Ly))


def polygon_sample(n, h, rho=None, R_circ=None, theta0=0.0, margin=3):
    """
    正 n 边形样品 (n>=3)。可用保角半径 rho 或外接圆半径 R_circ 指定大小(二选一)。
    推荐用 rho: 不同 n 的样品取相同 rho 时, 中心附近的“各向同性”热响应相同(见 docs/theory.md),
    只剩下 C_n 各向异性的差别, 便于公平比较。
    """
    if (rho is None) == (R_circ is None):
        raise ValueError("rho 与 R_circ 必须且只能给一个")
    if R_circ is None:
        R_circ = polygon_circumradius_from_rho(n, rho)
    else:
        rho = polygon_rho_from_circumradius(n, R_circ)
    half = R_circ + margin * h
    grid = Grid(half, half, h)
    inside = inside_regular_polygon(grid.X, grid.Y, n, R_circ, theta0)
    return _closed_shape_sample(grid, inside, f"poly{n}_rho{rho:g}",
                                dict(shape="polygon", n=n, rho=rho, R_circ=R_circ,
                                     theta0=theta0, h=h),
                                levelset_regular_polygon(n, R_circ, theta0))


def disk_sample(R, h, margin=3):
    """
    圆形样品, 半径 R, 边界恒温。作为“没有离散旋转对称性破缺”的对照组 (C_inf)。
    注意: 离散层面只剩网格的 C_4v 对称; 九点格式 + Gibou 边界使残余网格各向异性很小 (见 E0)。
    """
    grid = Grid(R + margin * h, R + margin * h, h)
    inside = inside_disk(grid.X, grid.Y, R)
    return _closed_shape_sample(grid, inside, f"disk_R{R:g}", dict(shape="disk", R=R, h=h),
                                levelset_disk(R))


def ellipse_sample(a, b, h, margin=3):
    grid = Grid(a + margin * h, b + margin * h, h)
    inside = inside_ellipse(grid.X, grid.Y, a, b)
    return _closed_shape_sample(grid, inside, f"ellipse_{a:g}x{b:g}",
                                dict(shape="ellipse", a=a, b=b, h=h), levelset_ellipse(a, b))


def surface_interface_sample(L=40.0, nd=-100.0, h=0.2, x_ext=2.0, margin=3):
    """
    导师 PRA 2008 "Nonlocal surface dipoles and vortices" 的几何:
      * 热介质:  -L/2 < x < 0,  |y| < L/2   (文中 eta -> 本代码 x, zeta -> y)
      * 线性介质: 0 < x < x_ext, 折射率差 n_d (<0), 不吸收、绝热 (界面 x=0 处 dn/dx = 0)
      * x = -L/2 及 y = +-L/2 三条边恒温 (n = 0); 光场在整个计算区外为 0。
    格心网格的节点在 x = +-h/2, ..., 界面 x=0 恰好位于两列节点正中, 零热流条件精确落在 x=0。
    线性介质只取到 x_ext (n_d=-100 时光场在其中的衰减长度约 0.07, x_ext=2 已足够)。
    """
    grid = Grid(L / 2 + margin * h, L / 2 + margin * h, h)
    X, Y = grid.X, grid.Y
    therm = (X > -L / 2) & (X < 0) & (np.abs(Y) < L / 2)
    linear = (X > 0) & (X < x_ext) & (np.abs(Y) < L / 2)
    th_type = np.full(grid.shape, SINK, dtype=np.int8)
    th_type[therm] = MEDIUM
    th_type[X > 0] = INSULATOR             # 界面右侧全部绝热(零热流)
    opt = therm | linear
    n_lin = np.where(linear, nd, 0.0)
    levelset = lambda X, Y: np.maximum(-L / 2 - X, np.abs(Y) - L / 2)   # 仅用于恒温的三条边
    return Sample(grid, th_type, opt, n_lin=n_lin, name=f"surface_L{L:g}_nd{nd:g}",
                  params=dict(shape="surface_interface", L=L, nd=nd, h=h, x_ext=x_ext),
                  levelset=levelset)


# ---------------------------------------------------------------------------
# 边界温度图样 (用于 E5: 热边界可重构调控)
# ---------------------------------------------------------------------------
def harmonic_multipole(n, Theta, R, theta0=0.0):
    """
    返回函数 f(X,Y) = Theta * (r/R)^n * cos(n*(theta - theta0))。
    它是 Laplace 方程的精确解; 在半径 R 的圆周上取值 Theta*cos(n(theta-theta0))。
    用作 SINK 点温度时, 样品内的“外加”折射率就是该多极调和场 (与功率无关的 C_n 各向异性)。
    n = 0 给出均匀温升 Theta (只平移传播常数)。
    """
    def f(X, Y):
        r = np.hypot(X, Y)
        th = np.arctan2(Y, X)
        return Theta * (r / R) ** n * np.cos(n * (th - theta0))
    f.description = f"Theta*(r/R)^n*cos(n(theta-theta0)), n={n}, Theta={Theta}, R={R}, theta0={theta0}"
    return f
