"""
operators.py —— 稀疏差分算子 (sparse finite-difference operators)
=================================================================

1) laplacian(mask, h, stencil)
   光场的 Laplace 算子, 掩膜外光场 = 0 (Dirichlet)。矩阵对称负定。

2) thermal_operator(th_type, kappa, th_bc, h, ...)
   热传导方程 div(kappa grad n) = -S 的离散:
       sum_{边 e} c_e * kappa_e * (n_nb - n_c) = -S_c
   邻居类型不同处理不同:
       MEDIUM    -> 正常耦合 (kappa_e 取两点热导率的调和平均)
       SINK      -> Dirichlet, 已知温度移到右端项 rhs (可用水平集做 Gibou 二阶修正, 见下)
       INSULATOR -> 该边热流为 0, 不出现 (Neumann)
   返回 (K, rhs), 离散方程写成  K @ n + rhs + S = 0。

模板 (stencil):
   stencil=5 : 标准五点格式, 截断误差含各向异性项 (h^2/12)(d^4/dx^4 + d^4/dy^4),
               即网格自带 C_4 对称的人为各向异性, 量级 (h k)^2。
   stencil=9 : 各向同性九点格式 (默认)
                   L9 = (2/3) L_+ + (1/3) L_x,
                   L_+ u = (u_E + u_W + u_N + u_S - 4u)/h^2,
                   L_x u = (u_NE + u_NW + u_SE + u_SW - 4u)/(2 h^2)   (对角方向, 间距 sqrt(2) h)
               其 h^2 阶误差是各向同性的 (h^2/12) Laplace^2, 网格各向异性被推到 O(h^4)。
               本课题研究的正是微弱的 C_n 各向异性 (n=3,4,5,6), 必须压低网格自身的 C_4 各向异性,
               否则网格效应可能冒充物理效应 (见 docs/numerics.md 中的检验)。

对称性很重要: 所有算子都按“边”组装, 每条边对两端点的贡献相同, 矩阵严格对称,
保证线性稳定性分析中的哈密顿结构 (本征值 ±lambda 成对), 不会出现离散化引入的虚假增长率。
"""

import numpy as np
import scipy.sparse as sp

from .geometry import MEDIUM, SINK

DEFAULT_STENCIL = 9


def _edges(stencil):
    """返回 [((dy, dx), 系数因子)], 边系数 = 因子 / h^2。"""
    nn = [((0, 1), 1.0), ((0, -1), 1.0), ((1, 0), 1.0), ((-1, 0), 1.0)]
    if stencil == 5:
        return nn
    if stencil == 9:
        return ([(d, 2.0 / 3.0) for d, _ in nn] +
                [((1, 1), 1.0 / 6.0), ((1, -1), 1.0 / 6.0), ((-1, 1), 1.0 / 6.0), ((-1, -1), 1.0 / 6.0)])
    raise ValueError("stencil 只能是 5 或 9")


def index_map(mask):
    """掩膜内的点编号 0..N-1 (按行优先, 即 C 顺序), 掩膜外为 -1。"""
    idx = -np.ones(mask.shape, dtype=np.int64)
    idx[mask] = np.arange(int(mask.sum()))
    return idx


def laplacian(mask, h, stencil=DEFAULT_STENCIL):
    """Laplace 算子 (Dirichlet: 掩膜外为 0)。返回 CSR 稀疏矩阵 (N x N)。"""
    idx = index_map(mask)
    N = int(mask.sum())
    iy, ix = np.nonzero(mask)
    Ny, Nx = mask.shape
    diag = np.zeros(N)
    rows, cols, vals = [], [], []
    for (dy, dx), f in _edges(stencil):
        c = f / h ** 2
        diag -= c
        jy, jx = iy + dy, ix + dx
        ok = (jy >= 0) & (jy < Ny) & (jx >= 0) & (jx < Nx)
        nb = np.full(N, -1, dtype=np.int64)
        nb[ok] = idx[jy[ok], jx[ok]]
        good = nb >= 0
        rows.append(np.arange(N)[good])
        cols.append(nb[good])
        vals.append(np.full(int(good.sum()), c))
    D = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N, N)) + sp.diags(diag)
    return D.tocsr()


def thermal_operator(th_type, kappa, th_bc, h, X=None, Y=None, levelset=None, bc_func=None,
                     stencil=DEFAULT_STENCIL, theta_min=0.01):
    """
    热传导算子。返回 (K, rhs), 方程为 K n + rhs + S = 0。

    若给出 levelset (及网格坐标 X, Y): 对每条“介质-热沉”边, 用水平集线性插值求真实边界在该边上的
    位置 (到介质点的距离 = theta * 边长, 0 < theta <= 1), 把该边的差商 (n_nb - n_c) 换成
    (n_b - n_c)/theta, n_b 为真实边界点上的温度 (bc_func 给出, 否则取 th_bc 在热沉点的值)。
    这就是 Gibou-Fedkiw-Cheng-Kang, J. Comput. Phys. 176, 205 (2002) 的对称二阶格式:
    只改对角元与右端, K 仍对称; 边界位置误差 O(h^2) (阶梯近似为 O(h), 会人为破坏 C_n 对称性)。
    theta 下限 theta_min 防止介质点紧贴边界时对角元过大。
    """
    med = th_type == MEDIUM
    idx = index_map(med)
    N = int(med.sum())
    iy, ix = np.nonzero(med)
    Ny, Nx = th_type.shape
    kc = kappa[iy, ix]
    diag = np.zeros(N)
    rhs = np.zeros(N)
    rows, cols, vals = [], [], []
    n_sink_edges = 0
    for (dy, dx), f in _edges(stencil):
        c = f / h ** 2
        jy, jx = iy + dy, ix + dx
        if ((jy < 0) | (jy >= Ny) | (jx < 0) | (jx >= Nx)).any():
            raise ValueError("热介质点紧贴网格边缘, 无法施加边界条件")
        t = th_type[jy, jx]
        kn = kappa[jy, jx]
        # --- 介质-介质: 调和平均热导率
        m = t == MEDIUM
        kf = 2.0 * kc[m] * kn[m] / (kc[m] + kn[m])
        rows.append(np.arange(N)[m])
        cols.append(idx[jy[m], jx[m]])
        vals.append(c * kf)
        diag[m] -= c * kf
        # --- 介质-热沉: Dirichlet, 已知值移到右端
        s = t == SINK
        theta = np.ones(int(s.sum()))
        Tb = th_bc[jy[s], jx[s]]
        if levelset is not None:
            xc, yc = X[iy[s], ix[s]], Y[iy[s], ix[s]]
            xn, yn = X[jy[s], jx[s]], Y[jy[s], jx[s]]
            pc, pn = levelset(xc, yc), levelset(xn, yn)
            with np.errstate(divide="ignore", invalid="ignore"):
                theta = np.where(pn > pc, pc / (pc - pn), 1.0)
            theta = np.clip(np.nan_to_num(theta, nan=1.0), theta_min, 1.0)
            if bc_func is not None:
                Tb = bc_func(xc + theta * (xn - xc), yc + theta * (yn - yc))
        diag[s] -= c * kc[s] / theta
        rhs[s] += c * kc[s] / theta * Tb
        n_sink_edges += int(s.sum())
        # --- 介质-绝热: 零热流, 什么都不加
    if n_sink_edges == 0:
        raise ValueError("热介质没有接触任何热沉(SINK), 纯 Neumann 问题奇异")
    K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N, N)) + sp.diags(diag)
    return K.tocsr(), rhs
