"""
E3c —— 圆盘中涡旋稳定性的径向一维独立复核（2026-10-13，第六轮文献调研之后）
==========================================================================
动机：Kartashov–Vysloukh–Torner 2007（Opt. Express 15, 9378）的正文声明圆柱热介质中 m ≤ 2 的涡旋在“全部存在域、任意 b”都稳定，
但该文没有给出 m=2 的稳定性图，也没有记录扫描范围（docs/literature_review.md Q1）。我们的二维计算（E3/E3b，方形网格
h=0.1/0.05）在低功率端给出 m=2 的窄失稳窗口（L=−2 分裂模，b R² < 31.8–32.2），m=1 也有窗口（L=−1）。
本脚本用与二维代码完全独立的方法复核：圆盘的转动对称性把问题化成径向一维——
  定态 q = f(r) e^{imθ} e^{ibξ}：牛顿迭代 + 沿 b 延拓；
  扰动 q = e^{ibξ}[w + U e^{λξ} + conj(V) e^{conj(λ)ξ}]，U = u(r) e^{i(m+L)θ}、V = v(r) e^{i(L−m)θ}、δN = ν(r) e^{iLθ} e^{λξ} + c.c.：
      iλ u = H_{m+L} u − f ν，  −iλ v = H_{L−m} v − f ν，  D_L ν = −f (u + v)，  H_k = −½ D_k − N + b，
  消去 ν（K = f G_L f，G_L = −D_L^{-1}，径向格林矩阵精确求逆）后 λ = −iμ，μ 为 [[H₁−K, −K], [K, −(H₂−K)]] 的本征值（稠密求解），
  增长率 Re λ = Im μ。通道 −L 的本征值是通道 L 的复共轭，所以只算 L = 0, −1, …, −(2m+1)。
没有笛卡尔网格、没有 C_4 网格误差、没有区域截断。径向格心网格 r_j = (j − ½)Δr、有限体积二阶格式（中心无通量、r=R 处 Dirichlet 镜像），
两档径向网格（M 与 2M）检验收敛。方程与归一化与 nlmodes 完全相同（i q_ξ = −½∇²q − N q，∇²N = −|q|²，N(R) = q(R) = 0）。

验证（云端原型，R=8）：L=0 通道的相位零模 |μ| = 8e-7（M=400）；m=2 的 U/b = 6.321（b=4，二维 6.31）；m=3 在 b=7 的增长率 0.2324（L=−3；二维 0.229）。

输出（results/E3_stability/<tag>/）：
  e3c_radial_growth.csv      每个 (m, b)：U、径向节点数（基本族为 0）、各通道 L 的最大增长率与对应的 Re μ（失稳模频率）、最不稳定通道
                             （各 m 先写分文件 e3c_radial_growth_m<m>.csv / e3c_radial_thresholds_m<m>.csv，再汇总；拆分运行用 --ms 与 --job）
  e3c_radial_thresholds.csv  每个 (m, L) 低功率窗口的上端 b_c（两档网格，二分到 --tol_b）、b_c R²（细网格）、两档是否一致（converged）
  e3c_radial_vs_2d.png       增长率 vs bR²：径向（线）与二维 E3/E3b（h=0.1、0.05，点）；虚线只画两档网格一致的端点（--plot_only 可单独重画）
运行：
  python e3c_radial_disk.py --quick     # 约 10 s（M=100/200，只检查流程）
  python e3c_radial_disk.py --hq        # 约 20–30 min（M=600/1200，论文数据；云端即可，不需要超算）
"""

import os
import sys
import time

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, pick, standard_args  # noqa: E402

from nlmodes import io_utils  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GROWTH_TOL = 1e-6        # 稠密求解时稳定态的 Im μ 是 ~1e-10 量级的舍入, 判据取 1e-6


def radial_grid(R, M):
    dr = R / M
    return (np.arange(M) + 0.5) * dr, dr


def radial_operator(r, dr, k):
    """有限体积 (1/r)(r f')' − k² f / r²：中心 r_{1/2} = 0 处无通量，外边界 f(R) = 0（镜像 f_{M+1} = −f_M）。稀疏三对角。"""
    rp, rm = r + dr / 2, r - dr / 2
    main = -(rp + rm) / (r * dr ** 2) - k ** 2 / r ** 2
    main[-1] -= rp[-1] / (r[-1] * dr ** 2)
    return sp.diags([rm[1:] / (r[1:] * dr ** 2), main, rp[:-1] / (r[:-1] * dr ** 2)], [-1, 0, 1], format="csc")


GUESSES = [(w, A) for w in (2.0, 3.0, 1.5, 1.0) for A in (2.0, 4.0, 1.0)]     # LG 型初值 A (r/w)^m exp(-r²/w²)


def radial_nodes(st):
    """f(r) 在 r < 0.95R 内的符号改变次数 (基本涡旋族为 0; 带径向节点的激发族 >= 1)。"""
    f = st["f"][st["r"] < 0.95 * st["r"][-1]]
    return int(np.sum(np.diff(np.sign(f[np.abs(f) > 1e-12 * np.abs(f).max()])) != 0))


def stationary(R, M, m, b, start=None, tol=1e-12, maxit=50):
    """½ D_m f + N f − b f = 0，D_0 N = −f²。start = 上一个解（延拓）或 None。
    start=None 时依次尝试 GUESSES 中的 LG 型初值, 只接受收敛、非零、无径向节点的解 (基本涡旋族): 通用初值在不同径向网格上
    可能收敛到零解或带节点的激发族 (第六轮: m=3 在 M=600 时落到 U(b=2)=58.9 的一支, 二维族为 20.6)。"""
    if start is None:
        r, dr = radial_grid(R, M)
        D0 = radial_operator(r, dr, 0)
        for w, A in GUESSES:
            f = A * (r / w) ** m * np.exp(-r ** 2 / w ** 2)
            st = stationary(R, M, m, b, start=dict(f=f, N=spla.spsolve(D0, -f ** 2)), tol=tol, maxit=maxit)
            if st["res"] < 1e-9 and st["U"] > 1e-3 and radial_nodes(st) == 0:
                return st
        raise RuntimeError(f"no nodeless m={m} state found at b={b} from the LG guesses")
    r, dr = radial_grid(R, M)
    Dm, D0 = radial_operator(r, dr, m), radial_operator(r, dr, 0)
    f, N = start["f"].copy(), start["N"].copy()
    res = np.inf
    for it in range(maxit):
        F1 = 0.5 * (Dm @ f) + N * f - b * f
        F2 = D0 @ N + f ** 2
        res = float(np.sqrt(np.sum(F1 ** 2 + F2 ** 2) * dr))
        if res < tol:
            break
        J = sp.bmat([[0.5 * Dm + sp.diags(N - b), sp.diags(f)], [sp.diags(2 * f), D0]], format="csc")
        d = spla.spsolve(J, -np.concatenate([F1, F2]))
        f, N = f + d[:M], N + d[M:]
    return dict(r=r, dr=dr, f=f, N=N, b=b, m=m, U=float(2 * np.pi * np.sum(r * f ** 2) * dr), res=res, its=it)


def channel_eigs(st, L):
    """通道 L 的 μ（λ = −iμ）。"""
    r, dr, f, N, b, m = st["r"], st["dr"], st["f"], st["N"], st["b"], st["m"]
    H1 = (-0.5 * radial_operator(r, dr, m + L) + sp.diags(b - N)).toarray()
    H2 = (-0.5 * radial_operator(r, dr, L - m) + sp.diags(b - N)).toarray()
    K = f[:, None] * (-la.inv(radial_operator(r, dr, L).toarray())) * f[None, :]
    return la.eigvals(np.block([[H1 - K, -K], [K, -(H2 - K)]]))


ZERO_MODE = 1e-4        # L=0 通道的相位零模是二重零本征值 (Jordan 块), 数值上分裂成 |μ| ~ 1e-5 的复数对, 不是失稳


def channel_growth(st, L):
    """(最大增长率, 该模的 |Re μ| = 失稳模频率)；稳定时频率记 nan。|μ| < ZERO_MODE 的近零模不计。"""
    mu = channel_eigs(st, L)
    mu = mu[np.abs(mu) >= ZERO_MODE]
    k = int(np.argmax(mu.imag))
    g = float(mu.imag[k])
    return (g, abs(float(mu.real[k]))) if g > GROWTH_TOL else (max(g, 0.0), float("nan"))


def regrid(st, R, M2):
    """把解插值到 M2 个格点的径向网格上 (作牛顿初值; 两档网格的收敛检验用)。"""
    r2, _ = radial_grid(R, M2)
    return dict(f=np.interp(r2, st["r"], st["f"]), N=np.interp(r2, st["r"], st["N"]), b=st["b"])


def continue_to(R, M, m, b, st):
    """从解 st 沿 b 以 ≤ 0.02 的步长延拓到 b。"""
    n = max(1, int(np.ceil(abs(b - st["b"]) / 0.02)))
    k0 = radial_nodes(st) if "r" in st else None
    for bb in np.linspace(st["b"], b, n + 1)[1:]:
        st = stationary(R, M, m, bb, start=st)
        if st["res"] > 1e-9:
            raise RuntimeError(f"Newton failed at m={m}, b={bb:.4f} (res {st['res']:.1e})")
        if k0 is not None and radial_nodes(st) != k0:
            raise RuntimeError(f"family switch at m={m}, b={bb:.4f}: radial nodes {k0} -> {radial_nodes(st)}")
    return st


def scan(R, M, m, bs):
    """从 b=2 的解出发，向下 / 向上延拓到 bs 中的每个点（升序返回）。"""
    st0 = stationary(R, M, m, 2.0)
    if st0["res"] > 1e-9:
        raise RuntimeError(f"no m={m} state at b=2 (res {st0['res']:.1e})")
    states = {}
    for side in (sorted(b for b in bs if b <= 2.0)[::-1], sorted(b for b in bs if b > 2.0)):
        st = st0
        for b in side:
            st = continue_to(R, M, m, b, st)
            states[b] = st
    return [states[b] for b in sorted(states)]


def window_edge(R, M, m, L, b_lo, b_hi, st_hi, tol_b):
    """b_lo 失稳、b_hi 稳定（通道 L）；二分到 tol_b，返回 (b_lo, b_hi)。"""
    st_lo = continue_to(R, M, m, b_lo, st_hi)
    assert channel_growth(st_lo, L)[0] > GROWTH_TOL and channel_growth(st_hi, L)[0] <= GROWTH_TOL
    while b_hi - b_lo > tol_b:
        bm = 0.5 * (b_lo + b_hi)
        st = continue_to(R, M, m, bm, st_hi)
        if channel_growth(st, L)[0] > GROWTH_TOL:
            b_lo = bm
        else:
            b_hi, st_hi = bm, st
    return b_lo, b_hi


def read_2d(tag, m):
    """二维 E3b（有通道标签）或 E3 的圆盘谱：[(b, growth)]；没有时返回 []。"""
    out = []
    for name, col in ((f"e3b_disk_m{m}_spectrum.csv", "growth"), (f"disk_m{m}_spectrum.csv", "max_growth")):
        f = os.path.join(io_utils.PROJECT_ROOT, "results", "E3_stability", tag, name)
        if os.path.exists(f):
            out = [(float(r["b"]), float(r[col])) for r in io_utils.read_csv(f)]
            break
    return out


def converged(t, rel=1e-3):
    """两档径向网格的窗口端点是否一致 (中点相对差 < rel)。"""
    a = 0.5 * (float(t["b_c_lo"]) + float(t["b_c_hi"]))
    b = 0.5 * (float(t["b_c2_lo"]) + float(t["b_c2_hi"]))
    return abs(a - b) < rel * b


def plot(out, R, ms, rows, thr):
    """增长率 vs bR²: 径向 (线, 每个失稳过的通道一条) 与二维 E3/E3b (点); 虚线 = 两档径向网格一致的窗口端点。
    只读 rows / thr (e3c_radial_growth.csv / e3c_radial_thresholds.csv 的行), --plot_only 可单独重画。"""
    fig, axs = plt.subplots(1, len(ms), figsize=(4.3 * len(ms), 3.4))
    axs = np.atleast_1d(axs)
    for ax, m in zip(axs, ms):
        mine = sorted((r for r in rows if int(float(r["m"])) == m), key=lambda r: float(r["b"]))
        bb = np.array([float(r["b"]) for r in mine])
        for L in range(0, -(2 * m + 2), -1):
            g = np.array([float(r[f"growth_L{L}"]) for r in mine])
            if (g > GROWTH_TOL).any():
                ax.plot(bb * R * R, g, "-", lw=1.4, label=f"radial, L={L}")
        if abs(R - 8.0) < 1e-9:
            for t, mk, lab in (("hq", "o", "2D, h=0.1"), ("hq_h005", "s", "2D, h=0.05")):
                d = read_2d(t, m)
                if d:
                    ax.plot([b * R * R for b, _ in d], [g for _, g in d], mk, ms=4.5, mfc="none", label=lab)
        for t in thr:
            if int(float(t["m"])) == m and converged(t):
                ax.axvline(float(t["bcR2"]), color="0.4", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("b R^2")
        ax.set_ylabel("growth rate Re(lambda)")
        ax.set_title(f"disk, m={m}: radial solver vs 2D grid", loc="left", fontsize=10)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "e3c_radial_vs_2d.png"), dpi=150)
    plt.close(fig)


def assemble(out, R, ms):
    """由各 m 的分文件 e3c_radial_growth_m<m>.csv / e3c_radial_thresholds_m<m>.csv 重建汇总表并作图 (拆分运行后用 --plot_only)。"""
    rows, thr = [], []
    for m in ms:
        fg, ft = (os.path.join(out, f"e3c_radial_{k}_m{m}.csv") for k in ("growth", "thresholds"))
        if os.path.exists(fg):
            rows += io_utils.read_csv(fg)
        if os.path.exists(ft):
            for t in io_utils.read_csv(ft):
                t["converged"] = converged(t)
                thr.append(t)
    have = sorted({int(float(r["m"])) for r in rows})
    io_utils.write_csv(os.path.join(out, "e3c_radial_growth.csv"), rows)
    io_utils.write_csv(os.path.join(out, "e3c_radial_thresholds.csv"), thr)
    if have:
        plot(out, R, have, rows, thr)
    for t in thr:
        print(f"m={t['m']} L={t['L']}: b_c R^2 = {float(t['bcR2']):.3f} (M={t['M2']}; M={t['M']}: "
              f"{0.5 * (float(t['b_c_lo']) + float(t['b_c_hi'])) * R * R:.3f}), max growth in the window "
              f"{float(t['max_growth_in_window']):.4f} at b={float(t['b_at_max']):g}"
              f"{'' if t['converged'] else '  [NOT converged between M and 2M]'}")
    return have


def main():
    p = standard_args(__doc__)
    p.add_argument("--R", type=float, default=8.0, help="圆盘半径（默认 8 = 二维图谱的保角半径 ρ）")
    p.add_argument("--ms", nargs="+", type=int, default=[1, 2, 3])
    p.add_argument("--M", type=int, default=None, help="径向格点数（默认 quick 100、main 400、hq 600；另用 2M 检验收敛）")
    p.add_argument("--tol_b", type=float, default=None, help="阈值二分的精度（默认 quick 2e-3、main/hq 1e-4）")
    p.add_argument("--plot_only", action="store_true",
                   help="不计算, 由各 m 的分文件重建 e3c_radial_growth.csv / _thresholds.csv 并重画对照图 (拆分运行后用)")
    args = finalize_args(p, "E3_stability")
    out = args.res_dir
    R = args.R
    if args.plot_only:
        print(f"assembled m = {assemble(out, R, args.ms)}; e3c_radial_vs_2d.png redrawn in {out}")
        return
    M = args.M or pick(args, quick=100, main=400, hq=600)
    tol_b = args.tol_b or pick(args, quick=2e-3, main=1e-4, hq=1e-4)
    db = pick(args, quick=0.1, main=0.05, hq=0.025)
    # 低功率端 (窗口所在) 用细步长 db, 之上逐渐放宽
    bs = sorted({round(b, 6) for b in np.concatenate([np.arange(db, 1.0, db), np.arange(1.0, 2.0 + 1e-9, 0.1),
                                                      np.arange(2.5, pick(args, quick=4.0, main=7.0, hq=7.0) + 1e-9, 0.5)])})
    print(f"radial check: R={R}, M={M} (dr={R / M:.4f}) and 2M, b = {bs[0]}..{bs[-1]} ({len(bs)} points), m = {args.ms}")
    t0 = time.time()
    rows, thr = [], []
    curves = {}
    for m in args.ms:
        Ls = list(range(0, -(2 * m + 2), -1))
        states = scan(R, M, m, bs)
        for st in states:
            gw = {L: channel_growth(st, L) for L in Ls}
            Lmax = max(Ls, key=lambda L: gw[L][0])
            row = dict(m=m, b=st["b"], bR2=st["b"] * R * R, U=st["U"], M=M, residual=st["res"], nodes=radial_nodes(st),
                       max_growth=gw[Lmax][0], L_max=Lmax if gw[Lmax][0] > GROWTH_TOL else "",
                       freq_max=gw[Lmax][1])
            for L in Ls:
                row[f"growth_L{L}"], row[f"freq_L{L}"] = gw[L]
            rows.append(row)
        curves[m] = (np.array([r["b"] for r in rows if r["m"] == m]),
                     {L: np.array([r[f"growth_L{L}"] for r in rows if r["m"] == m]) for L in Ls})
        io_utils.write_csv(os.path.join(out, f"e3c_radial_growth_m{m}.csv"), [r for r in rows if r["m"] == m])
        u2 = next((s["U"] for s in states if abs(s["b"] - 2.0) < 1e-9), float("nan"))
        print(f"  m={m}: {len(states)} states, U(b=2)={u2:.5f}, "
              f"max growth {max(r['max_growth'] for r in rows if r['m'] == m):.4f}  [{time.time() - t0:.0f}s]", flush=True)
        # 低功率窗口: 每个通道从低 b 起第一段失稳区间的上端 (之上至少有一个稳定的扫描点)
        bb, gl = curves[m]
        for L in Ls:
            un = gl[L] > GROWTH_TOL
            if not un[0] or un.all():
                continue
            k = int(np.argmin(un))                      # 第一个稳定点
            edges = []
            for MM in (M, 2 * M):
                st_hi = states[k] if MM == M else stationary(R, MM, m, bb[k], start=regrid(states[k], R, MM))
                if st_hi["res"] > 1e-9:
                    raise RuntimeError(f"Newton failed on the M={MM} grid at m={m}, b={bb[k]:.4f}")
                if channel_growth(st_hi, L)[0] > GROWTH_TOL:   # 细网格上该点仍失稳: 向上找稳定点
                    for b2 in bb[k + 1:]:
                        st_hi = continue_to(R, MM, m, b2, st_hi)
                        if channel_growth(st_hi, L)[0] <= GROWTH_TOL:
                            break
                lo, hi = window_edge(R, MM, m, L, bb[k - 1], st_hi["b"], st_hi, tol_b)
                edges.append((lo, hi))
                print(f"  m={m} L={L} M={MM}: window edge b_c in [{lo:.5f}, {hi:.5f}], b_c R^2 = {0.5 * (lo + hi) * R * R:.3f}"
                      f"  [{time.time() - t0:.0f}s]", flush=True)
            thr.append(dict(m=m, L=L, R=R, M=M, b_c_lo=edges[0][0], b_c_hi=edges[0][1], M2=2 * M,
                            b_c2_lo=edges[1][0], b_c2_hi=edges[1][1], bcR2=0.5 * sum(edges[1]) * R * R,
                            max_growth_in_window=float(gl[L][:k].max()),
                            b_at_max=float(bb[:k][np.argmax(gl[L][:k])])))
            thr[-1]["converged"] = converged(thr[-1])
        io_utils.write_csv(os.path.join(out, f"e3c_radial_thresholds_m{m}.csv"), [t for t in thr if t["m"] == m])

    # ---- 汇总 (含之前其他任务写的分文件) 与二维 E3/E3b 的对照图 (只有 R=8 时可比: 二维图谱的圆盘 ρ=8)
    assemble(out, R, [1, 2, 3, 4, 5, 6])
    print(f"E3c done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
