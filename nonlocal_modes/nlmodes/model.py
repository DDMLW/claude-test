"""
model.py —— 非局域热非线性模型 (把几何变成可计算的线性代数对象)
===================================================================

NonlocalModel(sample) 预先组装并分解好所有与解无关的矩阵:
    D      : 光场 Laplace 算子 (N_o x N_o), N_o = 光学未知点数
    K, rbc : 热传导算子与边界右端项 (N_t x N_t), N_t = 热介质点数
    Klu    : K 的稀疏 LU 分解 (每次求折射率只需回代, 很快)
    P      : 选择矩阵 (N_o x N_t): 热介质点上的折射率 -> 光学点
    n_lin  : 光学点上的线性折射率背景
    n_bc   : 仅由边界温度产生的折射率(调和部分, 与光强无关)

非线性折射率 (热介质点上):  n_t = n_bc + G[|q|^2],  G = (-K)^{-1} P^T  (正定, 非局域)
光所感受的总势:             V = n_lin + P n_t

stencil: 差分模板, 9 = 各向同性九点格式 (默认, 推荐), 5 = 五点格式 (仅用于对比/检验)。
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from .geometry import MEDIUM
from .operators import DEFAULT_STENCIL, index_map, laplacian, thermal_operator


class NonlocalModel:
    def __init__(self, sample, stencil=DEFAULT_STENCIL):
        self.sample = sample
        self.stencil = stencil
        g = sample.grid
        self.grid = g
        self.h = g.h
        self.opt_mask = sample.opt_mask
        self.th_mask = sample.th_type == MEDIUM
        self.No = int(self.opt_mask.sum())
        self.Nt = int(self.th_mask.sum())

        # 光学部分
        self.D = laplacian(self.opt_mask, self.h, stencil)
        self.n_lin = sample.n_lin[self.opt_mask].copy()
        self.xo = g.X[self.opt_mask]
        self.yo = g.Y[self.opt_mask]

        # 热学部分
        self.K, self.rbc = thermal_operator(sample.th_type, sample.kappa, sample.th_bc, self.h,
                                            X=g.X, Y=g.Y, levelset=sample.levelset,
                                            bc_func=sample.bc_func, stencil=stencil)
        self.Klu = splu(self.K.tocsc())

        # 光学点 <-> 热介质点 的对应关系
        idx_t = index_map(self.th_mask)
        oy, ox = np.nonzero(self.opt_mask)
        self.t_of_o = idx_t[oy, ox]                  # 光学点对应的热点编号, 非热介质为 -1
        sel = self.t_of_o >= 0
        self.P = sp.csr_matrix((np.ones(int(sel.sum())), (np.nonzero(sel)[0], self.t_of_o[sel])),
                               shape=(self.No, self.Nt))
        self.PT = self.P.T.tocsr()

        # 边界温度产生的(与光无关的)折射率
        self.n_bc = self.Klu.solve(-self.rbc) if np.any(self.rbc) else np.zeros(self.Nt)

    # ------------------------------------------------------------------
    # 折射率与势
    # ------------------------------------------------------------------
    def thermal_index(self, intensity):
        """intensity: 光学点上的 |q|^2 (长度 N_o) -> 热介质点上的折射率 n_t (含边界温度部分)。"""
        return self.n_bc + self.Klu.solve(-(self.PT @ intensity))

    def potential(self, n_t):
        """热介质点折射率 -> 光学点上的总势 V = n_lin + n。"""
        return self.n_lin + self.P @ n_t

    def potential_of_field(self, w):
        return self.potential(self.thermal_index(np.abs(w) ** 2))

    # ------------------------------------------------------------------
    # 常用量
    # ------------------------------------------------------------------
    def power(self, w):
        """能流 U = ∫∫ |q|^2 dx dy。"""
        return float(self.h ** 2 * np.sum(np.abs(w) ** 2))

    def to_grid(self, v, where="opt", fill=0.0):
        """把向量放回 (Ny,Nx) 网格。where='opt' 用光学掩膜, 'th' 用热介质掩膜。"""
        mask = self.opt_mask if where == "opt" else self.th_mask
        out = np.full(self.grid.shape, fill, dtype=np.result_type(v, float))
        out[mask] = v
        return out

    def from_grid(self, F, where="opt"):
        mask = self.opt_mask if where == "opt" else self.th_mask
        return np.asarray(F)[mask]

    def index_grid(self, n_t, fill=np.nan):
        """热介质点折射率 -> 网格 (非热介质点填 fill, 便于画图)。"""
        return self.to_grid(n_t, where="th", fill=fill)

    def residual(self, w, b):
        """定态方程残差 1/2 D w + (V - b) w (光学点上, 复数)。"""
        V = self.potential_of_field(w)
        return 0.5 * (self.D @ w) + (V - b) * w

    def describe(self):
        s = self.sample
        return (f"{s.name}: grid {s.grid.shape[1]}x{s.grid.shape[0]}, h={self.h:g}, stencil={self.stencil}, "
                f"optical unknowns={self.No}, thermal unknowns={self.Nt}")
