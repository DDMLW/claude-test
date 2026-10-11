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


def _augmented(model, sol):
    """增广线性化系统 A x = lambda B x 的稀疏矩阵 A 与 B 的对角线 (x = [a_r, a_i, dn])。"""
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
    return A, Bdiag, No, Nt


def eigenpairs_near(model, sol, sigma, k=12, tol=1e-10):
    """
    shift sigma 附近的 k 个本征值及其 Bogoliubov 模。返回按 |lambda - sigma| 排序的列表 [(lam, dict(u, v)), ...]
    (u, v 的定义见 stability_arnoldi; 用 mode_angular_content 可把每个本征值归入角向通道 L,
    用 krein_sign 得到 Krein 符号)。用于“按角向通道分解的 Bogoliubov 谱”(E3b)。
    """
    A, Bdiag, No, Nt = _augmented(model, sol)
    lu = splu((A - sigma * sp.diags(Bdiag)).astype(complex).tocsc())

    def mv(u):
        z = np.zeros(2 * No + Nt, complex)
        z[:2 * No] = u
        return lu.solve(z)[:2 * No]

    T = LinearOperator((2 * No, 2 * No), matvec=mv, dtype=complex)
    try:
        mu, vec = eigs(T, k=k, which="LM", return_eigenvectors=True, tol=tol)
    except ArpackNoConvergence as err:
        mu, vec = err.eigenvalues, err.eigenvectors
        if mu is None or len(mu) == 0:
            return []
    out = []
    for j in np.argsort(-np.abs(mu)):
        x = vec[:, j]
        p, s = x[:No], x[No:]
        u, v = 0.5 * (p + 1j * s), 0.5 * (p - 1j * s)
        scale = max(np.abs(u).max(), np.abs(v).max())
        out.append((complex(sigma + 1.0 / mu[j]), dict(u=u / scale, v=v / scale)))
    return out


def krein_sign(mode):
    """Krein 符号 sign(||u||^2 - ||v||^2): 纯虚本征值 (稳定振荡) 时有意义; 符号相反的两支相遇才会产生不稳定。"""
    d = np.vdot(mode["u"], mode["u"]).real - np.vdot(mode["v"], mode["v"]).real
    return int(np.sign(d))


def stability_arnoldi(model, sol, shifts=(0.05,), k=20, tol=1e-10, return_mode=False):
    """
    稀疏 shift-invert Arnoldi: 对每个 shift sigma 求离 sigma 最近的 k 个本征值。
    增广系统 A x = lambda B x,  x = [a_r, a_i, dn],  B = diag(I, I, 0):
        A = [[0, -L, -W_i P], [L, 0, W_r P], [2 P^T W_r, 2 P^T W_i, K]]
    算子 T u = [(A - sigma B)^{-1} [u; 0]]_{1:2N_o} 的本征值 mu = 1/(lambda - sigma)。
    返回所有找到的本征值 (去重前合并) 与 max Re(lambda)。
    return_mode=True: 另外返回增长率最大的本征模 mode = dict(lam, u, v) (Bogoliubov 形式,
        微扰 delta q = u exp(lambda xi) + conj(v) exp(conj(lambda) xi); 对电荷 m 的涡旋, u ~ e^{i(m+L)theta},
        v ~ e^{i(L-m)theta}, L 为微扰的角向指数, 见 mode_angular_content)。
    """
    t0 = time.time()
    A, Bdiag, No, Nt = _augmented(model, sol)
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
        # 释放本 shift 的 LU 再分解下一个: 否则新旧两份同时占内存 (h=0.05 时每份约 3.3 GB, 实测峰值 9.9 GB)
        del T, mv, lu
    lam = np.concatenate(found) if found else np.array([np.nan + 0j])
    out = dict(eigenvalues=lam, max_growth=float(np.nanmax(lam.real)) if found else float("nan"),
               n_failed_shifts=n_failed, time=time.time() - t0)
    if return_mode and found:
        out["mode"] = _eigenmode(A, Bdiag, No, Nt, lam[int(np.nanargmax(lam.real))], tol)
    return out


def _eigenmode(A, Bdiag, No, Nt, lam0, tol):
    """lam0 附近的本征矢 (shift 取 lam0 旁 1e-7 处, 避免矩阵奇异)。返回 dict(lam, u, v)。"""
    sigma = lam0 + 1e-7 * (1 + 1j)
    lu = splu((A - sigma * sp.diags(Bdiag)).astype(complex).tocsc())

    def mv(u):
        z = np.zeros(2 * No + Nt, complex)
        z[:2 * No] = u
        return lu.solve(z)[:2 * No]

    T = LinearOperator((2 * No, 2 * No), matvec=mv, dtype=complex)
    try:
        mu, vec = eigs(T, k=1, which="LM", return_eigenvectors=True, tol=tol)
    except ArpackNoConvergence as err:
        if err.eigenvalues is None or len(err.eigenvalues) == 0:
            return None
        mu, vec = err.eigenvalues[:1], err.eigenvectors[:, :1]
    x = vec[:, 0]
    p, s = x[:No], x[No:]
    u, v = 0.5 * (p + 1j * s), 0.5 * (p - 1j * s)
    scale = max(np.abs(u).max(), np.abs(v).max())
    return dict(lam=complex(sigma + 1.0 / mu[0]), u=u / scale, v=v / scale)


def mode_angular_content(model, mode, m, center=(0.0, 0.0), r=None, nsym=None, lmax=None, w=None):
    """
    Bogoliubov 模 (u, v) 在半径 r 圆周上的角向谐波分析。微扰 delta q = u e^{lambda xi} + conj(v) e^{conj(lambda) xi};
    对电荷 m 的涡旋 (圆盘中) u ~ e^{i(m+L)theta}, v ~ e^{i(L-m)theta} (L 为微扰角向指数);
    多边形中各自还混有 l + j n 的谐波。r 默认取定态 w 的亮环半径 (不要取模式自身的峰: 例如分裂模 u ~ e^{i0theta}
    在中心最大, 而 v ~ r^{|L-m|} 在中心几乎为 0, 在小半径处会读错 v 的谐波)。返回 dict:
        l_u, l_v      u、v 的主导角向指数;  w_u, w_v  其权重 (占该圆周上总功率的比例)
        L             由 u 得到的微扰角向指数 l_u - m;  L_from_v = l_v + m (两者相等说明结构自洽)
        kappa         L 对 n 取模 (落在 (-n/2, n/2]); kappa = 0 表示微扰保持 C_n 对称 (只改变同一不可约表示内的结构)
        norm_ratio    ||v|| / ||u||;  top_u, top_v  前三个谐波 [(l, 权重), ...]
    例: 圆盘中 m=2 的低功率失稳为 L = -2 (u ~ e^{0 i theta}, v ~ e^{-4 i theta}) —— 双荷涡旋分裂为两个单位涡旋的四极模。
    """
    from . import diagnostics as dg
    ug = dg.field_grid(model, mode["u"])
    vg = dg.field_grid(model, mode["v"])
    if r is None:
        ref = dg.field_grid(model, w) if w is not None else np.abs(ug) + np.abs(vg) + 0j
        r = dg.ring_radius(ref, model.grid, center)
    lmax = lmax or (abs(m) + 12)
    ls, cu = dg.angular_spectrum(ug, model.grid, center, r, lmax=lmax)
    _, cv = dg.angular_spectrum(vg, model.grid, center, r, lmax=lmax)
    wu = np.abs(cu) ** 2 / max(np.sum(np.abs(cu) ** 2), 1e-300)
    wv = np.abs(cv) ** 2 / max(np.sum(np.abs(cv) ** 2), 1e-300)
    iu, iv = int(np.argmax(wu)), int(np.argmax(wv))
    L = int(ls[iu]) - m
    out = dict(r=float(r), l_u=int(ls[iu]), l_v=int(ls[iv]), w_u=float(wu[iu]), w_v=float(wv[iv]), L=L,
               L_from_v=int(ls[iv]) + m,
               norm_ratio=float(np.linalg.norm(mode["v"]) / max(np.linalg.norm(mode["u"]), 1e-300)),
               top_u=[(int(ls[j]), round(float(wu[j]), 3)) for j in np.argsort(wu)[::-1][:3]],
               top_v=[(int(ls[j]), round(float(wv[j]), 3)) for j in np.argsort(wv)[::-1][:3]])
    if nsym:
        k = L % nsym
        out["kappa"] = int(k - nsym if k > nsym / 2 else k)
    return out
