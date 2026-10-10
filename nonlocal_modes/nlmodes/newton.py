"""
newton.py —— 定态解的牛顿迭代 (Newton iteration for stationary nonlinear modes)
=============================================================================

求  q(x,y,xi) = w(x,y) exp(i b xi)  形式的定态解, w 一般为复数 (涡旋) 或实数 (偶极/多极):

    F1 = 1/2 D w_r + (V - b) w_r = 0
    F2 = 1/2 D w_i + (V - b) w_i = 0          V = n_lin + P n_t
    F3 = K n_t + rbc + P^T (w_r^2 + w_i^2) = 0

与导师 MATLAB 代码 (Relax2Dlocal.m) 的关系:
  * 那份代码是“局域”非线性 (折射率只依赖当地光强), 雅可比矩阵只含对角非线性项;
  * 这里折射率由热传导方程决定 (非局域), 我们把 n_t 也当作未知量, 雅可比矩阵仍然是稀疏的
    3x3 分块矩阵, 用稀疏 LU 直接求解, 不需要任何稠密矩阵:

        J = [[ L ,  0 , W_r P ],
             [ 0 ,  L , W_i P ],
             [2P^T W_r, 2P^T W_i, K]] ,    L = 1/2 D + diag(V - b)

  * 每步迭代后用 n_t = G[|w|^2] 精确更新折射率 (F3 恒为 0), 等价于对“消去 n 后的非局域方程”做
    精确牛顿迭代, 收敛性与 MATLAB 代码一样是二次收敛。
  * 复数解有 U(1) 相位不变性 (w -> w e^{i phi}), 雅可比矩阵奇异。这里把 |w| 最大点 p 的虚部
    固定为 0 ("钉住相位", 替换方程 F2[p])。因为离散 Laplace 算子对称, 被替换掉的方程在收敛时
    自动满足 (Noether 恒等式 sum(w_r F2 - w_i F1) = 0)。

可选“定功率”模式 (fixed_power=U0): 把 b 也当未知量, 增加约束 U(w) = U0, 用加边 (bordering)
算法求解, 只需同一个 LU 分解做两次回代。适合按功率扫描、或比较不同尺寸的样品。
"""

import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu


def _is_real(w):
    w = np.asarray(w)
    return np.max(np.abs(w.imag)) <= 1e-13 * max(1.0, np.max(np.abs(w.real)))


def _replace_row_with_unit(J, row):
    """把稀疏矩阵 J (COO) 的第 row 行换成单位行 e_row (用于钉住相位)。"""
    J = J.tocoo()
    keep = J.row != row
    rows = np.concatenate([J.row[keep], [row]])
    cols = np.concatenate([J.col[keep], [row]])
    vals = np.concatenate([J.data[keep], [1.0]])
    return sp.csc_matrix((vals, (rows, cols)), shape=J.shape)


def jacobian(model, w, V, b, real):
    """组装牛顿雅可比矩阵 (见文件头部公式)。"""
    L = 0.5 * model.D + sp.diags(V - b)
    wr = w.real
    WrP = sp.diags(wr) @ model.P
    PTWr2 = 2.0 * (model.PT @ sp.diags(wr))
    if real:
        return sp.bmat([[L, WrP], [PTWr2, model.K]], format="coo")
    wi = w.imag
    WiP = sp.diags(wi) @ model.P
    PTWi2 = 2.0 * (model.PT @ sp.diags(wi))
    return sp.bmat([[L, None, WrP], [None, L, WiP], [PTWr2, PTWi2, model.K]], format="coo")


def align_phase(w, ref):
    """把 w 乘上一个全局相位, 使其与 ref 的重叠积分为正实数 (连续追踪解族时用)。"""
    s = np.vdot(w, ref)
    if abs(s) == 0:
        return w
    return w * (s / abs(s))


GLOBALIZATIONS = ("armijo", "none", "nonmonotone")

# 默认收敛判据 (max|F|)。实验脚本按精度档位修改它: main 1e-9, hq 1e-12 (见 experiments/common.py)。
# 不能小于 ~1e-13: 双精度机器精度约 2.2e-16, 乘上差分算子系数 (~1/h^2) 后残差的舍入噪声即为这个量级。
DEFAULT_TOL = 1e-9


def solve(model, w0, b, real=None, fixed_power=None, tol=None, maxit=40, min_power=1e-6,
          verbose=False, globalization="armijo"):
    """
    牛顿迭代求定态解。

    参数
    ----
    model       : NonlocalModel
    w0          : 初始猜测 (光学点上的向量, 实数或复数)
    b           : 传播常数 (定功率模式下作为 b 的初值)
    real        : True 只解实数解 (偶极/多极); None 时根据 w0 自动判断
    fixed_power : None 或目标能流 U0 (此时 b 也是未知量)
    tol         : 收敛判据 max|F| < tol (定功率时同时要求 |U-U0|/U0 < tol); None 时用 DEFAULT_TOL。
                  若残差已降到 1e-10 以下且牛顿步不再使其明显下降 (已到舍入误差底线), 也判为收敛
                  (只在 hq 档 tol=1e-12 时起作用)。
    globalization: 'armijo' (单调回溯线搜索, 默认) | 'none' (纯牛顿, 步长恒为 1) |
                   'nonmonotone' (与最近 5 步价值函数的最大值比较, 允许暂时上升)。
                   激发态(涡旋)的牛顿收敛域很小, 冷启动时三种策略各有成败, 见 solve_robust。

    返回 dict: w, n_t, b, power, converged, iterations, res_hist, message, time
    """
    t0 = time.time()
    if tol is None:
        tol = DEFAULT_TOL
    w = np.asarray(w0, dtype=complex).copy()
    if real is None:
        real = _is_real(w)
    No = model.No
    pin = None
    if real:
        w = w.real.astype(complex)
    else:
        pin = int(np.argmax(np.abs(w)))
        w *= np.exp(-1j * np.angle(w[pin]))     # 使 w[pin] 为正实数
    b = float(b)
    hist = []
    message = "max iterations reached"
    converged = False

    def evaluate(w, b):
        """返回 n_t, V, 残差 R, 收敛判据 res (最大模), 线搜索用的价值函数 merit (2-范数)。
        牛顿方向一定是 2-范数平方的下降方向, 但不一定是最大模的下降方向, 所以线搜索用 merit。"""
        n_t = model.thermal_index(np.abs(w) ** 2)
        V = model.potential(n_t)
        R = 0.5 * (model.D @ w) + (V - b) * w
        res = float(np.max(np.abs(R)))
        # 线搜索价值函数: 去掉被“钉住”的那个方程 Im F[pin] (钉相位后的牛顿步只对其余方程是精确牛顿步;
        # 由 Noether 恒等式, 其余方程满足时 Im F[pin] 自动为 0)
        Rm = R
        if pin is not None:
            Rm = R.copy()
            Rm[pin] = R.real[pin]
        merit = float(np.linalg.norm(Rm)) * model.h
        if fixed_power is not None:
            rp = abs(model.power(w) - fixed_power) / fixed_power
            res = max(res, rp)
            merit = float(np.hypot(merit, rp))
        return n_t, V, R, res, merit

    n_t, V, R, res, merit = evaluate(w, b)
    merit_hist = [merit]
    it = 0
    for it in range(maxit + 1):
        hist.append(res)
        P = model.power(w)
        if verbose:
            print(f"    newton it={it:2d}  res={res:.3e}  b={b:.6f}  U={P:.6f}")
        if not np.isfinite(res):
            message = "NaN/Inf encountered"
            break
        if P < min_power:
            message = "converged to (or near) the zero solution"
            break
        if res < tol:
            converged = True
            message = "converged"
            break
        if it >= 3 and res < 1e-10 and res > 0.5 * hist[-2]:
            converged = True
            message = f"converged (round-off floor reached, residual {res:.1e})"
            break
        if it == maxit:
            break
        if it >= 15 and res > 1e-2 * hist[0]:
            message = "stagnation (residual not reduced 100x after 15 iterations)"
            break

        # ---- 组装并分解雅可比矩阵
        J = jacobian(model, w, V, b, real)
        if real:
            F = np.concatenate([R.real, np.zeros(model.Nt)])
        else:
            J = _replace_row_with_unit(J, No + pin)
            F = np.concatenate([R.real, R.imag, np.zeros(model.Nt)])
            F[No + pin] = w.imag[pin]
        try:
            lu = splu(J.tocsc(), permc_spec="COLAMD")
        except RuntimeError as err:          # 奇异矩阵
            message = f"singular Jacobian: {err}"
            break
        y1 = -lu.solve(F)
        db = 0.0
        if fixed_power is not None:
            # 加边算法: J dx + c db = -F,  d^T dx = -(U - U0)
            c = np.zeros_like(F)
            c[:No] = -w.real
            if not real:
                c[No:2 * No] = -w.imag
                c[No + pin] = 0.0
            y2 = lu.solve(c)
            d = np.zeros_like(F)
            d[:No] = 2 * model.h ** 2 * w.real
            if not real:
                d[No:2 * No] = 2 * model.h ** 2 * w.imag
            FP = model.power(w) - fixed_power
            db = (d @ y1 + FP) / (d @ y2)
            y1 = y1 - db * y2
        dw = y1[:No].astype(complex)
        if not real:
            dw = dw + 1j * y1[No:2 * No]

        # ---- 回溯线搜索 (Armijo 条件: 价值函数必须下降, 否则步长减半; 实在不降则判定失败)
        ref = max(merit_hist[-5:]) if globalization == "nonmonotone" else merit
        alpha = 1.0
        while True:
            w_try = w + alpha * dw
            b_try = b + alpha * db
            n_try, V_try, R_try, res_try, merit_try = evaluate(w_try, b_try)
            if globalization == "none" or merit_try < (1.0 - 1e-4 * alpha) * ref:
                break
            alpha *= 0.5
            if alpha < 1.0 / 256:
                break
        if alpha < 1.0 / 256:
            message = "line search failed (no descent direction; try a better initial guess)"
            break
        w, b, n_t, V, R, res, merit = w_try, b_try, n_try, V_try, R_try, res_try, merit_try
        merit_hist.append(merit)

    return dict(w=w, n_t=n_t, b=float(b), power=model.power(w), converged=converged,
                iterations=it, res_hist=hist, message=message, real=bool(real),
                time=time.time() - t0, globalization=globalization)


def solve_robust(model, w0, b, strategies=GLOBALIZATIONS, **kw):
    """依次尝试几种全局化策略 (默认 armijo -> none -> nonmonotone), 返回第一个收敛的结果 (或最后一次的失败结果)。"""
    sol = None
    for g in strategies:
        sol = solve(model, w0, b, globalization=g, **kw)
        if sol["converged"]:
            break
    return sol


def continuation(model, w_start, p_start, p_stop, dp, param="b", real=None, tol=None,
                 maxit=20, max_halvings=3, max_rel_jump=0.5, callback=None, verbose=True,
                 strategies=("armijo", "none")):
    """
    自然参数延拓 (natural-parameter continuation), 沿 b 或 U(=power) 追踪一个解族。

    param      : 'b' -> 固定 b 求解; 'P' -> 固定能流 U 求解 (b 为未知量)
    p_start    : 起始参数值; p_stop: 终止值; dp: 步长 (带符号, 例如向下扫 b 用负数)
    max_rel_jump: 新解与预测值的相对偏差上限, 超过则认为“跳到了别的解支”, 步长减半重试
    callback(sol): 每得到一个收敛解调用一次 (用于保存数据/画图)

    返回 (sols, info): sols 为按顺序的解列表; info 记录停止原因 (例如“存在阈值”)。
    """
    def run(w_guess, pval, b_guess):
        # 有好的预测初值时牛顿通常 3-6 步收敛; 20 步还不收敛基本说明走到了解族端点 (阈值/分岔)
        if param == "b":
            return solve_robust(model, w_guess, pval, strategies=strategies, real=real, tol=tol,
                                maxit=maxit)
        return solve_robust(model, w_guess, b_guess, strategies=strategies, real=real,
                            fixed_power=pval, tol=tol, maxit=maxit)

    sols = []
    first = run(w_start, p_start, p_start if param == "b" else _rayleigh_b(model, w_start))
    if not first["converged"]:
        return sols, dict(stop_reason=f"initial solve failed: {first['message']}",
                          last_param=p_start)
    first["param"] = p_start
    sols.append(first)
    if callback:
        callback(first)
    if verbose:
        print(f"  [{param}={p_start:.4f}] b={first['b']:.5f} U={first['power']:.5f} "
              f"its={first['iterations']}")

    p = p_start
    step = dp
    direction = np.sign(dp)
    halvings = 0
    stop_reason = "reached end of range"

    def next_grid_point(pv):
        # 参数网格 p_start + k*dp 上沿延拓方向的下一个格点。步长减半后的点落在格点之间,
        # 恢复步长时先回到格点上, 使解族的参数值与是否发生过减半无关 (不同网格的结果可逐点比较)
        k = np.floor((pv - p_start) / dp + 1e-9) + 1
        return p_start + k * dp

    while direction * (p_stop - p) > 1e-12:
        p_new = p + step
        g = next_grid_point(p)
        if direction * (p_new - g) > -1e-12:
            p_new = g
        if direction * (p_new - p_stop) > 0:
            p_new = p_stop
        cur = sols[-1]
        if len(sols) >= 2:
            prev = sols[-2]
            w_prev = align_phase(prev["w"], cur["w"])
            s = (p_new - cur["param"]) / (cur["param"] - prev["param"])
            w_guess = cur["w"] + s * (cur["w"] - w_prev)
            b_guess = cur["b"] + s * (cur["b"] - prev["b"])
        else:
            w_guess, b_guess = cur["w"], cur["b"]
        sol = run(w_guess, p_new, b_guess)
        ok = sol["converged"]
        if ok:
            wa = align_phase(sol["w"], w_guess)
            jump = np.linalg.norm(wa - w_guess) / max(np.linalg.norm(w_guess), 1e-30)
            if jump > max_rel_jump:
                ok = False
                sol["message"] = f"branch jump suspected (rel. change {jump:.2f})"
        if ok:
            sol["param"] = p_new
            sols.append(sol)
            if callback:
                callback(sol)
            if verbose:
                print(f"  [{param}={p_new:.4f}] b={sol['b']:.5f} U={sol['power']:.5f} "
                      f"its={sol['iterations']}")
            p = p_new
            halvings = 0
            if abs(step) < abs(dp):
                step = np.sign(step) * min(abs(dp), 2 * abs(step))
        else:
            halvings += 1
            if verbose:
                print(f"  [{param}={p_new:.4f}] failed ({sol['message']}); halving step")
            if halvings > max_halvings:
                stop_reason = (f"continuation stopped near {param}={p:.5f}: "
                               f"{sol['message']} (possible existence threshold/bifurcation)")
                break
            step = step / 2
    return sols, dict(stop_reason=stop_reason, last_param=p)


def _rayleigh_b(model, w):
    """由 Rayleigh 商估计 b = <w,(1/2 D + V) w>/<w,w> (定功率模式的初值)。"""
    V = model.potential_of_field(w)
    num = np.vdot(w, 0.5 * (model.D @ w) + V * w).real
    return float(num / np.vdot(w, w).real)
