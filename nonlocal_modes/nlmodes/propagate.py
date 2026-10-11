"""
propagate.py —— 传播动力学 (beam propagation)
=============================================

    i dq/dxi = -1/2 D q - V(q) q,     V = n_lin + P n_t,   K n_t + rbc + P^T |q|^2 = 0

数值格式: Strang 分步 (二阶精度)
    q <- exp(i dxi/2 V) q               (非线性半步, 不改变 |q|, 因此 V 不变)
    q <- (I - i dxi/4 D)^{-1} (I + i dxi/4 D) q     (线性 Crank-Nicolson 整步, 酉算子, 严格保能流)
    重新解热传导方程得到新的 V
    q <- exp(i dxi/2 V) q
每步只需: 1 次复稀疏 LU 回代 + 1 次热传导回代。能流守恒到机器精度 (可作为程序正确性检查)。

重要: 空间离散 (D, K) 与 newton.py 完全一致, 因此牛顿解在此格式下是严格的定态
(只有 O(dxi^2) 的时间离散误差), 稳定性测试不会因为“离散不匹配”而出现虚假的振荡。

边界温度随传播距离变化 (E5: 热边界图样的切换/旋转) 通过 nbc_func(xi) 传入,
它返回热介质点上“仅由边界温度产生的折射率” (见 model.n_bc)。
"""

import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from . import diagnostics as dg


def add_noise(w, variance=0.01, seed=0):
    """输入噪声: q(0) = w (1 + rho), rho 为复高斯白噪声, <|rho|^2> = variance (论文常用 0.01)。"""
    rng = np.random.default_rng(seed)
    s = np.sqrt(variance / 2.0)
    rho = s * (rng.standard_normal(w.shape) + 1j * rng.standard_normal(w.shape))
    return w * (1.0 + rho)


def propagate(model, q0, z_max, dz=0.02, diag_every=50, snapshot_zs=(), nbc_func=None,
              w_ref=None, nsym=None, center=(0.0, 0.0), track_singularities=False,
              sing_region_threshold=1e-3, sing_core_radius=None, verbose=True, progress_every=None,
              custom_diag=None, sing_every=1):
    """
    传播 q0 (光学点上的复向量) 到 xi = z_max。

    参数
    ----
    dz            : 步长 (默认 0.02; 必须做步长收敛检验)
    diag_every    : 每隔多少步记录一次诊断量
    snapshot_zs   : 需要保存完整场的传播距离列表
    nbc_func      : None 或 函数 xi -> 热介质点上的边界温度折射率 (默认用 model.n_bc, 不随 xi 变)
    w_ref         : 参考定态解 (用于计算偏离度 deviation = || |q|-|w_ref| || / ||w_ref||)
    nsym          : 若给出, 记录 C_nsym 不可约表示权重 (角赝动量守恒检验)
    track_singularities: 记录光束区域内的单位奇点 (中心半径 sing_core_radius 以外, 默认 1.25h),
                  以及中心荷 central (以 center 为圆心、半径 0.5h 的圆周绕数) —— 动态嬗变的核心观测量
    sing_every    : 奇点检测较慢, 每隔几次诊断才做一次 (其余诊断时刻记为 None)
    custom_diag   : None 或 函数 (z, q_grid) -> dict, 在每次诊断时调用, 结果按时间顺序存入 out["custom"]

    返回 dict:
        z, power, peak, deviation, xc, yc, (irreps), (singularities), snapshots{z: q}, time
    """
    t0 = time.time()
    No = model.No
    q = np.asarray(q0, complex).copy()
    Id = sp.identity(No, format="csc", dtype=complex)
    Am = (Id - 0.25j * dz * model.D).tocsc()
    Ap = (Id + 0.25j * dz * model.D).tocsr()
    lu = splu(Am)

    def potential(q, z):
        n_t = model.Klu.solve(-(model.PT @ (np.abs(q) ** 2)))
        n_t += model.n_bc if nbc_func is None else nbc_func(z)
        return model.n_lin + model.P @ n_t

    nsteps = int(round(z_max / dz))
    snap_steps = {int(round(zs / dz)): zs for zs in snapshot_zs}
    rec = dict(z=[], power=[], peak=[], deviation=[], xc=[], yc=[])
    if nsym:
        rec["irreps"] = []
    if track_singularities:
        rec["singularities"] = []
        rec["central"] = []
        core = 1.25 * model.h if sing_core_radius is None else sing_core_radius
    if custom_diag is not None:
        rec["custom"] = []
    n_diag = [0]
    snapshots = {}
    ref_norm = None if w_ref is None else np.linalg.norm(w_ref)

    def record(z, q):
        rec["z"].append(z)
        rec["power"].append(model.power(q))
        rec["peak"].append(float(np.abs(q).max()))
        xc, yc = dg.center_of_mass(model, q)
        rec["xc"].append(xc)
        rec["yc"].append(yc)
        if w_ref is not None:
            rec["deviation"].append(float(np.linalg.norm(np.abs(q) - np.abs(w_ref)) / ref_norm))
        if nsym or track_singularities or custom_diag is not None:
            qg = dg.field_grid(model, q)
            if nsym:
                rec["irreps"].append(dg.cn_weights(qg, model.grid, nsym, center))
            if track_singularities:
                if n_diag[0] % max(sing_every, 1) == 0:
                    reg = dg.beam_region(qg, model.grid, rel_threshold=sing_region_threshold) & model.opt_mask
                    rec["singularities"].append(dg.find_singularities(qg, model.grid, reg, exclude_center=center,
                                                                      exclude_radius=core))
                else:
                    rec["singularities"].append(None)
                rec["central"].append(dg.winding_on_circle(qg, model.grid, center, 0.5 * model.h))
            if custom_diag is not None:
                rec["custom"].append(custom_diag(z, qg))
        n_diag[0] += 1

    V = potential(q, 0.0)
    record(0.0, q)
    if 0 in snap_steps:
        snapshots[0.0] = q.copy()
    if progress_every is None:
        progress_every = max(1, nsteps // 10)
    for step in range(1, nsteps + 1):
        z = step * dz
        q *= np.exp(0.5j * dz * V)
        q = lu.solve(Ap @ q)
        V = potential(q, z)
        q *= np.exp(0.5j * dz * V)
        if step % diag_every == 0 or step == nsteps:
            record(z, q)
        if step in snap_steps:
            snapshots[snap_steps[step]] = q.copy()
        if verbose and step % progress_every == 0:
            msg = f"    xi={z:9.2f}  U={rec['power'][-1]:.8f}  peak={rec['peak'][-1]:.5f}"
            if w_ref is not None:
                msg += f"  dev={rec['deviation'][-1]:.3e}"
            print(msg + f"  ({time.time() - t0:.0f}s)")
    out = {k: (np.array(v) if k not in ("irreps", "singularities", "custom") else v) for k, v in rec.items()}
    out["snapshots"] = snapshots
    out["final"] = q
    out["time"] = time.time() - t0
    return out
