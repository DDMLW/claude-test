"""
stability.py —— 线性稳定性分析 (linear stability analysis)
==========================================================

微扰 q = [w + a(x,y) e^{lambda xi}] e^{i b xi}, a = a_r + i a_i, 线性化得到

    lambda a_r = -L a_i - w_i dn
    lambda a_i =  L a_r + w_r dn,            L = 1/2 D + diag(V - b)
    dn = G[ 2 (w_r a_r + w_i a_i) ],         G = (-K)^{-1} (非局域热响应, 对称正定)

即 lambda x = J H x (哈密顿结构): 本征值以 (lambda, -lambda, lambda*, -lambda*) 成组出现。
存在 Re(lambda) > 0 即不稳定, 增长率 = max Re(lambda)。

提供两种方法:
  1) stability_dense  : 在光束附近的圆形区域(半径 crop_radius)内构造稠密矩阵, 求全部本征值。
                        可靠 (不会漏掉本征值), 适合 N_s <= ~6000 个格点。
                        截断区域外微扰取 0; 非局域响应 dn 仍按整个样品计算 (G 由全样品热算子给出)。
                        结论必须做 crop_radius 收敛检验 (见 E3 README)。
                        实数解 (偶极/多极) 自动用 lambda^2 = eig(-L L_+) 降维加速。
  2) stability_arnoldi: 整个样品上的稀疏增广系统 + shift-invert Arnoldi, 求离若干个 shift 最近的本征值。
                        适合大网格的“定点精算”(先用方法 1 在粗网格上确定不稳定本征值的位置)。
"""

import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigs, splu

from .diagnostics import center_of_mass


def _crop_indices(model, center, crop_radius):
    if crop_radius is None:
        return np.arange(model.No)
    r = np.hypot(model.xo - center[0], model.yo - center[1])
    return np.nonzero(r < crop_radius)[0]


def green_block(model, S, chunk=400):
    """G_S = (-K)^{-1} 在光学点集合 S 上的稠密块 (S 中非热介质点对应的行列为 0)。"""
    Ns = len(S)
    t = model.t_of_o[S]
    th = np.nonzero(t >= 0)[0]
    G = np.zeros((Ns, Ns))
    for start in range(0, len(th), chunk):
        c = th[start:start + chunk]
        rhs = np.zeros((model.Nt, len(c)))
        rhs[t[c], np.arange(len(c))] = 1.0
        X = model.Klu.solve(-rhs)               # (-K)^{-1} e_j
        G[np.ix_(th, c)] = X[t[th], :]
    return G


def stability_dense(model, sol, crop_radius=None, center=None, max_points=6000, return_modes=False):
    """
    稠密本征值法。返回 dict:
        eigenvalues : 全部本征值 (复数数组)
        max_growth  : max Re(lambda)
        Ns, crop_radius, center, time
        (return_modes=True 时另含 mode: 增长率最大的本征矢对应的微扰场 a (光学点向量, 截断区外为 0))
    """
    t0 = time.time()
    w = sol["w"]
    b = sol["b"]
    if center is None:
        center = center_of_mass(model, w)
    S = _crop_indices(model, center, crop_radius)
    Ns = len(S)
    if Ns > max_points:
        raise ValueError(f"截断区域内有 {Ns} 个格点 > max_points={max_points}; 请减小 crop_radius 或放大 h")
    V = model.potential(model.thermal_index(np.abs(w) ** 2))
    Lfull = (0.5 * model.D + sp.diags(V - b)).tocsr()
    L = Lfull[S][:, S].toarray()
    G = green_block(model, S)
    wr, wi = w.real[S], w.imag[S]
    real = np.max(np.abs(wi)) <= 1e-12 * np.max(np.abs(wr))
    mode = None
    if real:
        Lp = L + 2.0 * wr[:, None] * G * wr[None, :]
        A = -L @ Lp                                   # lambda^2 的本征值问题
        if return_modes:
            import scipy.linalg as sla
            mu, vec = sla.eig(A)
        else:
            mu = np.linalg.eigvals(A)
        lam = np.sqrt(mu.astype(complex))
        lam = np.where(lam.real < 0, -lam, lam)         # 取 Re >= 0 的那支
        eig_all = np.concatenate([lam, -lam])
        if return_modes:
            k = int(np.argmax(lam.real))
            ar = vec[:, k]
            ai = (Lp @ ar) / lam[k] if abs(lam[k]) > 0 else np.zeros_like(ar)   # lambda a_i = L_+ a_r
            mode = np.zeros(model.No, complex)
            mode[S] = ar + 1j * ai
    else:
        GWr = G * wr[None, :]
        GWi = G * wi[None, :]
        A11 = -2.0 * wi[:, None] * GWr
        A12 = -L - 2.0 * wi[:, None] * GWi
        A21 = L + 2.0 * wr[:, None] * GWr
        A22 = 2.0 * wr[:, None] * GWi
        M = np.block([[A11, A12], [A21, A22]])
        if return_modes:
            import scipy.linalg as sla
            eig_all, vec = sla.eig(M)
            k = int(np.argmax(eig_all.real))
            mode = np.zeros(model.No, complex)
            mode[S] = vec[:Ns, k] + 1j * vec[Ns:, k]
        else:
            eig_all = np.linalg.eigvals(M)
    out = dict(eigenvalues=eig_all, max_growth=float(np.max(eig_all.real)), Ns=int(Ns),
               crop_radius=crop_radius, center=tuple(center), time=time.time() - t0, real=bool(real))
    if return_modes:
        out["mode"] = mode
    return out


def stability_arnoldi(model, sol, shifts=(0.05,), k=20, tol=1e-10):
    """
    稀疏 shift-invert Arnoldi: 对每个 shift sigma 求离 sigma 最近的 k 个本征值。
    增广系统 A x = lambda B x,  x = [a_r, a_i, dn],  B = diag(I, I, 0):
        A = [[0, -L, -W_i P], [L, 0, W_r P], [2 P^T W_r, 2 P^T W_i, K]]
    算子 T u = [(A - sigma B)^{-1} [u; 0]]_{1:2N_o} 的本征值 mu = 1/(lambda - sigma)。
    返回所有找到的本征值 (去重前合并) 与 max Re(lambda)。
    """
    t0 = time.time()
    w = sol["w"]
    b = sol["b"]
    No, Nt = model.No, model.Nt
    V = model.potential(model.thermal_index(np.abs(w) ** 2))
    L = 0.5 * model.D + sp.diags(V - b)
    Wr, Wi = sp.diags(w.real), sp.diags(w.imag)
    A = sp.bmat([[None, -L, -Wi @ model.P],
                 [L, None, Wr @ model.P],
                 [2 * model.PT @ Wr, 2 * model.PT @ Wi, model.K]], format="csc")
    Bdiag = np.concatenate([np.ones(2 * No), np.zeros(Nt)])
    found = []
    n_failed = 0          # ARPACK 未完全收敛的 shift 个数 (只用已收敛的本征值, 不中断长任务)
    for sigma in shifts:
        try:
            lu = splu((A - sigma * sp.diags(Bdiag)).astype(complex).tocsc())
        except RuntimeError:                 # shift 恰好是本征值 (矩阵奇异): 稍微移开
            sigma = sigma + 1e-6 * (1 + 1j)
            lu = splu((A - sigma * sp.diags(Bdiag)).astype(complex).tocsc())

        def mv(u, lu=lu):
            z = np.zeros(2 * No + Nt, complex)
            z[:2 * No] = u
            return lu.solve(z)[:2 * No]

        T = LinearOperator((2 * No, 2 * No), matvec=mv, dtype=complex)
        try:
            mu = eigs(T, k=k, which="LM", return_eigenvectors=False, tol=tol)
        except ArpackNoConvergence:
            try:                             # 加大 Krylov 子空间重试一次
                mu = eigs(T, k=k, which="LM", return_eigenvectors=False, tol=tol,
                          ncv=min(2 * No - 1, max(4 * k + 1, 60)))
            except ArpackNoConvergence as err:
                mu = err.eigenvalues
                n_failed += 1
        if len(mu):
            found.append(sigma + 1.0 / mu)
    lam = np.concatenate(found) if found else np.array([np.nan + 0j])
    return dict(eigenvalues=lam, max_growth=float(np.nanmax(lam.real)) if found else float("nan"),
                n_failed_shifts=n_failed, time=time.time() - t0)
