"""
E5d —— 加热态 (默认四极加热双涡旋) 在固定能流处的低频 Bogoliubov 谱 vs Theta
==========================================================================
目的 (第四轮, theory.md §5.1): “受迫分裂模”图像中, 旋转图样以 Omega 转动时在 2 Omega = omega_{-2}(Theta)
(L=-2、Krein + 的分裂模 = 双涡旋整体转动的模) 处共振, 预言同向旋转的临界角速度 Omega_c ≈ omega_{-2}/2;
二阶 (∝ Theta^2) 驱动作用在 L=-4 通道上, 共振在 4 Omega = omega_{-4}。本脚本给出这些频率随 Theta 的变化,
与 E5b 的 Omega 扫描 (hq_omega) 对照; Theta 增大到 U 处失稳时也能看到哪一支模先变软 / 碰撞。

做法: 半径 R 的圆盘, 边界温度 Theta cos(n theta), 拓扑荷 m, 固定能流 U (与 E5b 相同); 从 Theta=0 沿 Theta 延拓
(割线预测, 细步长); 每个 Theta 在 Im(lambda) = 0.02 ... --wmax 之间放 shift, 求附近的本征值, 去掉 |lambda| 很小的零模,
按角向通道 L (u ~ e^{i(m+L)theta}) 与 Krein 符号分类 (加热态中 L 只是主导谐波, kappa = L mod n 才是好量子数)。

输出 (results/E5_thermal_reconfiguration/<tag>/):
    e5d_n<n>_m<m>_U<U>_modes.csv   每个 (Theta, 本征值): b、lam_re、lam_im、L、L_from_v、kappa、krein、w_u
    e5d_n<n>_m<m>_U<U>_modes.png   频率 vs Theta (颜色 = L, 实心 = Krein +, 空心 = Krein -; 红框 = 不稳定)
    e5d_n<n>_m<m>_U<U>_summary.csv 每个 Theta: b、omega_{-2} (转动/分裂模, 从 Theta=0 的 L=-2 Krein + 最低支连续追踪) 与
                                   一阶共振 omega_{-2}/2、L=-4 模与二阶共振 omega_{-4}/4、最软的破坏 C_n 的 Krein - 模
                                   omega_odd (实本征值时记 0)、最大增长率、由 omega_odd^2 线性外推的 Theta_c
                                   (--table_only 由 _modes.csv 重建; 2026-10-13 起按连续性追踪, 不再按 L 标签)
运行:
  python e5d_twin_modes.py --quick                       (约 1 分钟)
  python e5d_twin_modes.py --hq --job e5d_n2              (单核约 0.5-1 h)
  python e5d_twin_modes.py --hq --fields_only --thetas 0 0.25 0.5 1 1.2 --job e5d_fields
                                                         (只导出定态场 results/<tag>/fields/, 供 Fig.6 插图; 几分钟)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, heated_disk, initial_solution, param_continuation, pick, standard_args  # noqa: E402

from nlmodes import guesses, io_utils, newton, stability  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# 通道的固定配色 (与 analysis/style.py 同一色盲安全调色板, 固定顺序; |L| > 4 归入灰色)
L_COLORS = {-4: "#2a78d6", -3: "#eb6834", -2: "#1baf7a", -1: "#eda100", 0: "#e87ba4", 1: "#008300", 2: "#6250d6",
            3: "#e34948"}


def track(rows, thetas, start_L, krein, kappa0=True, w_jump=0.05):
    """
    沿 Theta 连续追踪一支模: Theta 最小处取 L=start_L、给定 Krein 符号的最低支, 之后每个 Theta 取同一 Krein 符号
    (kappa0=True 时只看 kappa=0) 中频率最接近上一点的模 (相差 < w_jump)。不能只按 L 标签选: 强加热时主导谐波会变
    (第四轮: Theta>=1 时转动/分裂模的 u 主谐波从 e^{0i theta} 变成 e^{-2i theta}, 标签由 L=-2 变成 -4, 按标签会跳到另一支 0.825)。
    返回 {Theta: 频率} (没有找到的 Theta 不在其中)。
    """
    out, prev = {}, None
    for T in thetas:
        cand = [r for r in rows if r["Theta"] == T and r["krein"] == krein and r["lam_re"] < 1e-5
                and (not kappa0 or r["kappa"] == 0)]
        if prev is None:
            cand = sorted((r["lam_im"] for r in cand if r["L"] == start_L))
            if cand:
                out[T] = prev = cand[0]
            continue
        if cand:
            best = min(cand, key=lambda r: abs(r["lam_im"] - prev))
            if abs(best["lam_im"] - prev) < w_jump:
                out[T] = prev = best["lam_im"]
    return out


def summarize(rows, thetas, n):
    """每个 Theta: 转动/分裂模 omega_{-2} 与 L=-4 模 omega_{-4} (连续追踪)、最软的破坏 C_n 的 Krein - 模 omega_odd
    (实本征值时记 0)、最大增长率; Theta_c 由稳定点 omega_odd^2 对 Theta 的线性外推给出。"""
    w2 = track(rows, thetas, -2, 1)
    w4 = track(rows, thetas, -4, 1)
    summ = []
    for T in thetas:
        mine = [r for r in rows if r["Theta"] == T]
        odd = sorted(r["lam_im"] for r in mine if r["kappa"] != 0 and r["krein"] == -1 and r["lam_re"] < 1e-5)
        gmax = max([r["lam_re"] for r in mine] + [0.0])
        summ.append(dict(Theta=T, b=mine[0]["b"] if mine else "", omega_m2=w2.get(T, ""),
                         Omega_c_pred=w2[T] / 2 if T in w2 else "", omega_m4=w4.get(T, ""),
                         Omega_res2=w4[T] / 4 if T in w4 else "",
                         omega_odd=(0.0 if gmax > 1e-5 else (odd[0] if odd else "")), max_growth=gmax,
                         n_modes=len(mine)))
    st = [(r["Theta"], r["omega_odd"] ** 2) for r in summ if r["omega_odd"] not in ("", 0.0) and r["max_growth"] < 1e-5]
    if len(st) >= 3:
        fit = np.polyfit([t for t, _ in st], [w for _, w in st], 1)
        for r in summ:
            r["Theta_c_fit"] = -fit[1] / fit[0] if fit[0] < 0 else ""
    return summ


def main():
    p = standard_args(__doc__)
    p.add_argument("--R", type=float, default=10.0)
    p.add_argument("--n", type=int, default=2, help="加热图样阶数 (默认 2: 四极加热双涡旋)")
    p.add_argument("--m", type=int, default=2)
    p.add_argument("--U", type=float, default=15.0, help="固定能流 (与 E5b 相同)")
    p.add_argument("--thetas", nargs="+", type=float, default=None, help="默认 hq: 0-1.5 共 10 个")
    p.add_argument("--wmax", type=float, default=0.62, help="shift 的最高频率 (Im)")
    p.add_argument("--table_only", action="store_true", help="不计算, 由 _modes.csv 重建 _summary.csv (连续追踪)")
    p.add_argument("--fields_only", action="store_true",
                   help="只求各 Theta 的定态并导出裁剪后的场 (fields/<名>_T<Theta>.npz, 单精度, 供论文插图), 不算谱")
    args = finalize_args(p, "E5_thermal_reconfiguration")
    out = args.res_dir
    R, h, n, m, U = args.R, args.h, args.n, args.m, args.U
    if args.table_only:
        name = f"e5d_n{n}_m{m}_U{U:g}"
        rows = [{k: (float(v) if k not in ("L", "L_from_v", "kappa", "krein") else int(float(v))) for k, v in r.items()}
                for r in io_utils.read_csv(os.path.join(out, f"{name}_modes.csv"))]
        thetas = sorted({r["Theta"] for r in rows})
        io_utils.write_csv(os.path.join(out, f"{name}_summary.csv"), summarize(rows, thetas, n))
        print(f"{name}_summary.csv rebuilt from {len(rows)} modes ({len(thetas)} Theta)")
        return
    thetas = sorted(args.thetas if args.thetas is not None else
                    pick(args, quick=[0.0, 0.5], main=[0.0, 0.25, 0.5, 0.75, 1.0],
                         hq=[0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.1, 1.2, 1.3, 1.5]))
    shifts = [1j * w for w in np.arange(0.02, args.wmax + 1e-9, pick(args, quick=0.15, main=0.08, hq=0.06))]
    kk = pick(args, quick=6, main=8, hq=10)
    name = f"e5d_n{n}_m{m}_U{U:g}"

    # ---- Theta=0 的定态 (固定能流), 再沿 Theta 延拓 (细步长 + 割线预测, 与 E5b 相同)
    M0 = heated_disk(R, h, n, 0.0)
    first = initial_solution(M0, lambda W: guesses.lg_vortex(M0, m, W), [3.0, 3.5, 2.5], verbose=True)
    s0 = newton.solve_robust(M0, first["w"], first["b"], fixed_power=U)
    targets = sorted({round(x, 6) for T in thetas if T > 0 for x in np.linspace(0.0, T, 9)[1:]})
    scan = param_continuation(lambda v: heated_disk(R, h, n, v), 0.0, M0, s0, targets,
                              lambda model, w0, prev: newton.solve_robust(model, w0, prev["b"], fixed_power=U, maxit=25),
                              secant=True, verbose=False)
    states = {0.0: (M0, s0)}
    for T, Mt, st in scan:
        states[round(T, 6)] = (Mt, st)
    if args.fields_only:
        fdir = io_utils.ensure_dir(os.path.join(out, "fields"))
        for T in thetas:
            if round(T, 6) not in states:
                continue
            Mt, st = states[round(T, 6)]
            q = Mt.to_grid(np.asarray(st["w"], complex))
            g = Mt.grid
            half = 4.0                                    # 双涡旋的亮环半径约 2.3: |x|,|y| <= 4 包含整个光束
            ix, iy = np.where(np.abs(g.x) <= half)[0], np.where(np.abs(g.y) <= half)[0]
            sl = (slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1))
            fn = os.path.join(fdir, f"{name}_T{T:g}.npz")
            np.savez_compressed(fn, x=g.x[sl[1]].astype(np.float32), y=g.y[sl[0]].astype(np.float32),
                                psi=q[sl].astype(np.complex64), b=st["b"], U=U, h=g.h, Theta=T, n=n, m=m)
            print(f"Theta={T:g}: b={st['b']:.4f}, field -> {fn} ({os.path.getsize(fn) / 1024:.0f} kB)", flush=True)
        print("\nE5d done (fields only).")
        return

    rows = []
    for T in thetas:
        if round(T, 6) not in states:
            print(f"Theta={T:g}: not reached by the continuation (stopped at {scan[-1][0] if scan else 0:g})",
                  flush=True)
            continue
        Mt, st = states[round(T, 6)]
        seen = []
        for sg in shifts:
            for lam, md in stability.eigenpairs_near(Mt, st, sg, k=kk):
                if lam.imag < -1e-9 or abs(lam) < 1e-5 or any(abs(lam - x) < 1e-6 for x in seen):
                    continue
                seen.append(lam)
                ac = stability.mode_angular_content(Mt, md, m, w=st["w"])
                rows.append(dict(Theta=T, b=st["b"], lam_re=lam.real, lam_im=lam.imag, L=ac["L"],
                                 L_from_v=ac["L_from_v"], kappa=((ac["L"] + n // 2) % n) - n // 2,
                                 krein=stability.krein_sign(md), w_u=round(ac["w_u"], 3)))
        summ = summarize(rows, sorted({r["Theta"] for r in rows}), n)
        cur = summ[-1]
        print(f"Theta={T:g}: b={st['b']:.4f}, {cur['n_modes']} modes, omega_rot={cur['omega_m2']}, "
              f"omega(L=-4,K+)={cur['omega_m4']}, omega_odd={cur['omega_odd']}, max growth {cur['max_growth']:.2e}",
              flush=True)
        io_utils.write_csv(os.path.join(out, f"{name}_modes.csv"), rows)
        io_utils.write_csv(os.path.join(out, f"{name}_summary.csv"), summ)

    # ---- 图: 频率 vs Theta
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for L in sorted({r["L"] for r in rows}):
        col = L_COLORS.get(L, "0.6")
        for kr, mk in ((1, "o"), (-1, "^")):
            sel = [r for r in rows if r["L"] == L and r["krein"] == kr]
            if sel:
                ax.plot([r["Theta"] for r in sel], [r["lam_im"] for r in sel], mk, ms=5, ls="none", color=col,
                        mfc=col if kr == 1 else "none", mec=col,
                        label=f"L={L:+d}" if kr == 1 or not [r for r in rows if r["L"] == L and r["krein"] == 1]
                        else None)
    uns = [r for r in rows if r["lam_re"] > 1e-5]
    ax.plot([r["Theta"] for r in uns], [r["lam_im"] for r in uns], "s", ms=10, mfc="none", mec="#e34948", ls="none",
            label="unstable")
    ax.set_xlabel("heating amplitude Theta")
    ax.set_ylabel("Im(lambda)")
    ax.set_title(f"{name}: low-lying Bogoliubov modes (o Krein +, ^ Krein -)", fontsize=8)
    ax.legend(fontsize=6, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{name}_modes.png"), dpi=130)
    plt.close(fig)
    print("\nE5d done.")


if __name__ == "__main__":
    main()
