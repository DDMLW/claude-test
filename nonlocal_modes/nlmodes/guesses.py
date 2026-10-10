"""
guesses.py —— 牛顿迭代的初始猜测 (initial guesses)
===================================================

牛顿法只在初值足够好时收敛, 激发态(涡旋、多极)尤其如此。这里提供:
  * lg_vortex      : 拉盖尔-高斯型涡旋 (r/w0)^|m| exp(-r^2/w0^2) exp(i m theta)
  * multipole_ring : 环上 k 个高斯光斑, 相邻光斑相位差可设 (偶极、四极、项链等)
  * fit_amplitude  : 给定形状, 用 Rayleigh 商匹配传播常数 b 求振幅
  * best_guess     : 在一组宽度 w0 中挑选残差最小者 (简单的“变分”初值)
"""

import numpy as np


def lg_vortex(model, m, w0, center=(0.0, 0.0), amp=1.0, theta0=0.0):
    """光学点上的 LG 型涡旋初值, 拓扑荷 m (可为负; m=0 即高斯基模)。"""
    x = model.xo - center[0]
    y = model.yo - center[1]
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    return amp * (r / w0) ** abs(m) * np.exp(-(r / w0) ** 2) * np.exp(1j * m * (th - theta0))


def multipole_ring(model, k, ring_radius, spot_width, phase_step=np.pi, center=(0.0, 0.0),
                   theta0=0.0, amp=1.0):
    """
    环上 k 个高斯光斑: 第 j 个位于角度 theta0 + 2 pi j/k, 相位 j*phase_step。
    k=2, phase_step=pi -> 偶极;  k=4, phase_step=pi -> 四极;  phase_step=2 pi m/k -> 离散涡旋。
    """
    x = model.xo - center[0]
    y = model.yo - center[1]
    out = np.zeros(model.No, complex)
    for j in range(k):
        a = theta0 + 2 * np.pi * j / k
        xc, yc = ring_radius * np.cos(a), ring_radius * np.sin(a)
        out += np.exp(1j * j * phase_step) * np.exp(-((x - xc) ** 2 + (y - yc) ** 2) / spot_width ** 2)
    out *= amp
    if np.allclose(out.imag, 0, atol=1e-12 * np.abs(out).max()):
        return out.real.astype(complex)
    return out


def fit_amplitude(model, shape, b):
    """
    对给定形状 phi 求振幅 A, 使 Rayleigh 商等于 b:
        <phi,(1/2 D + n_lin + n_bc) phi> + A^2 <phi, n_1 phi> = b <phi,phi>,  n_1 = G[|phi|^2]。
    返回 A (若无正实解返回 None, 说明这个 b 对该形状太小/太大)。
    """
    phi = shape
    n1 = model.P @ (model.thermal_index(np.abs(phi) ** 2) - model.n_bc)   # 单位振幅光强产生的折射率
    V0 = model.potential(model.n_bc)                                     # 与光无关的势 (线性背景 + 边界温度)
    num = b * np.vdot(phi, phi).real - np.vdot(phi, 0.5 * (model.D @ phi) + V0 * phi).real
    den = np.vdot(phi, n1 * phi).real
    if den <= 0 or num <= 0:
        return None
    return float(np.sqrt(num / den))


def best_guess(model, make_shape, b, widths):
    """
    make_shape(w0) -> 光学点上的形状 (单位振幅)。对每个 w0 用 fit_amplitude 定振幅,
    计算定态方程残差 ||F||, 返回残差最小的初值 (以及对应 w0)。
    """
    best = (np.inf, None, None)
    for w0 in widths:
        phi = make_shape(w0)
        A = fit_amplitude(model, phi, b)
        if A is None:
            continue
        w = A * phi
        res = np.linalg.norm(model.residual(w, b)) / np.linalg.norm(w)
        if res < best[0]:
            best = (res, w, w0)
    if best[1] is None:
        raise ValueError("所有宽度下都找不到正的振幅, 请调整 b 或宽度范围")
    return best[1], best[2]


def guess_at_power(model, make_shape, power, widths):
    """给定能流 U 的初值: 振幅由 U 固定, 在 widths 中挑残差(以 Rayleigh 商为 b)最小者。"""
    best = (np.inf, None, None, None)
    for w0 in widths:
        phi = make_shape(w0)
        w = phi * np.sqrt(power / model.power(phi))
        V = model.potential_of_field(w)
        b = np.vdot(w, 0.5 * (model.D @ w) + V * w).real / np.vdot(w, w).real
        res = np.linalg.norm(0.5 * (model.D @ w) + (V - b) * w) / np.linalg.norm(w)
        if res < best[0]:
            best = (res, w, w0, b)
    return best[1], best[2], best[3]


# ---------------------------------------------------------------------------
# 自洽线性模式初值 (适用于强非局域介质: 非线性模式 ≈ 自感应势中的线性模式)
# ---------------------------------------------------------------------------
def linear_modes(model, V, k=8):
    """
    势 V 中线性算子 1/2 D + V 的 k 个最高本征值 (即束缚得最紧的模式), 按本征值降序返回 (vals, vecs)。
    用 shift-invert (shift 取在全部本征值之上), 对几万个未知量只需 1-2 秒。
    """
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    H = (0.5 * model.D + sp.diags(V)).tocsc()
    vals, vecs = eigsh(H, k=k, sigma=float(V.max()) + 0.1, which="LM")
    order = np.argsort(-vals)
    return vals[order], vecs[:, order]


def scf_guess(model, power, combine, intensity0, iters=15, mix=0.5, k=8):
    """
    自洽迭代: 由当前光强算折射率 -> 求线性模式 -> combine(vecs, vals) 组合出新场 -> 与旧光强混合。
    power     : 固定能流 U
    combine   : 函数 (vecs, vals) -> 光学点上的复场, 例如 lambda v, e: v[:,1] + 1j*0.5*v[:,2] (非正则涡旋)
    intensity0: 初始光强 (光学点向量)
    返回 (w, b): 场与 Rayleigh 商估计的传播常数, 可直接作为定功率牛顿迭代的初值。
    """
    w = np.sqrt(np.clip(intensity0, 0, None)).astype(complex)
    w *= np.sqrt(power / model.power(w))
    for _ in range(iters):
        V = model.potential_of_field(w)
        vals, vecs = linear_modes(model, V, k=k)
        new = np.asarray(combine(vecs, vals), complex)
        new *= np.sqrt(power / model.power(new))
        w = np.sqrt((1 - mix) * np.abs(w) ** 2 + mix * np.abs(new) ** 2) * np.exp(1j * np.angle(new))
        w *= np.sqrt(power / model.power(w))
    V = model.potential_of_field(w)
    b = np.vdot(w, 0.5 * (model.D @ w) + V * w).real / np.vdot(w, w).real
    return w, float(b)


# ---------------------------------------------------------------------------
# 理论指导的“嬗变”初值: 按选择定则放置中心奇点与卫星环
# ---------------------------------------------------------------------------
def transmuted_guess(model, m, n, w0, r_s=1.0, theta_s=0.0, ellipticity=1.0, center=(0.0, 0.0)):
    """
    按 theory.selection_rule(m, n) 构造初值: 中心荷 v (z^v 或 conj(z)^|v|) 乘以各卫星环上的因子 (z - z_k),
    再乘高斯包络 exp(-(x^2/(e w0)^2 + y^2/w0^2))。第 k 个环 (含 N 个单位卫星, N 为 n 的整数倍) 的卫星放在
    半径 r_s*(1 + 0.6 k_sub) 处, 方位角 theta_s + 2 pi j / n。强各向异性 (例如矩形中的 m>=3) 时
    正则 LG 初值离真解太远, 这种初值往往能收敛。
    """
    from .theory import selection_rule
    rule = selection_rule(m, n)
    z = (model.xo - center[0]) + 1j * (model.yo - center[1])
    v = rule["central"]
    field = z ** v if v >= 0 else np.conj(z) ** (-v)
    field = field.astype(complex) / w0 ** abs(v)
    ring_index = 0
    for (count, charge, _) in rule["rings"]:
        for k_sub in range(count // n):
            rad = r_s * (1.0 + 0.6 * ring_index)
            ring_index += 1
            for j in range(n):
                a = theta_s + 2 * np.pi * j / n
                zk = rad * np.exp(1j * a)
                field *= ((z - zk) if charge > 0 else np.conj(z - zk)) / w0
    env = np.exp(-((model.xo - center[0]) / (ellipticity * w0)) ** 2 - ((model.yo - center[1]) / w0) ** 2)
    return field * env
