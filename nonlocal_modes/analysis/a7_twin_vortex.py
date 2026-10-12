"""
a7_twin_vortex.py —— 四极 (C_2) 加热下的 "双涡旋": 可调、稳定、可旋转、可释放 (论文 Fig.6)
==========================================================================================
只读取已上传的结果表 (圆盘 R=10, h=0.1, U=15):
  E5a  results/E5_thermal_reconfiguration/<tag>/heating_n2_m2.csv          卫星半径 r_s vs Theta
  E5c  results/E5_thermal_reconfiguration/<tag>/e5c_windows.csv, e5c_summary.csv   (2:2) 的失稳窗口与 U=15 处的增长率
  E4b  results/E4_dynamic_transmutation/<tag>/e4b_heat2_m2_U15*_timeseries.csv     LG m=2 射入加热样品后双涡旋的形成
  E5b  results/E5_thermal_reconfiguration/<tag>/control_n2_T{0.5,1}_*_tracks.csv   旋转跟随与关断后的自由共转
  E3b  results/E3_stability/hq_hib/e3b_disk_m2_channels.csv                       圆盘 L=-2 (Krein +) 通道频率 (预言共转角速度)

六幅子图:
  (a) 双涡旋半间距 r_s vs Theta (对数): 渐近 r_s ∝ Theta^(1/2) (牛顿多边形: n=2, m=2); 2h 为可分辨下限
  (b) (Theta, b) 平面的失稳窗口: 低 b 窗口的上端随 Theta 线性上移; U=15 的工作点在所有 Theta <= 1 都稳定,
      上端外推到 b(U=15) 给出这一功率下可调范围的上限 Theta_c
  (c) LG m=2 光束射入加热样品: 双涡旋在 xi≈4 内形成并保持 (C_2 对称性保持到 xi=600)
  (d) 加热图样以 Omega 旋转时, 双涡旋轴相对图样的滞后角: 慢转 (0.005) 绝热跟随; 快转 (0.05) Theta=0.5 跟随,
      Theta=1 在 xi≈164 失锁 (成核出涡旋-反涡旋对)
  (e) xi=150 关断加热后的自由共转: 角速度 = 圆盘 m=2 涡旋 L=-2 (Krein +) Bogoliubov 模频率的一半
      (u ∝ e^{i(m+L)theta} = 常数, 在中心处填充 -> 两个单位涡旋, 其轴以 omega/2 旋转)
  (f) 双涡旋间距 d(xi)
输出: results/analysis/<tag>/twin_vortex.png、twin_vortex_summary.csv
运行:  python a7_twin_vortex.py
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import style  # noqa: E402

ROOT = os.path.join(io_utils.PROJECT_ROOT, "results")
R_SAMPLE = 10.0          # E5 圆盘半径
RHO_CHANNELS = 8.0       # E3b 通道谱所用圆盘半径 (标度律 b -> b R^2, lambda -> lambda R^2 换算到 R=10)
H_GRID = 0.1
OMEGA = {"rotate_slow": 0.005, "rotate_fast": 0.05}   # e5b_dynamic_control.py 中的图样角速度 (hq)


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def e5_file(tag, name):
    return os.path.join(ROOT, "E5_thermal_reconfiguration", tag, name)


def pair_track(tag, theta, case):
    """两个 +1 奇点的轴角 (mod pi) 与间距; 奇点数不等于 2 (或出现 -1 荷) 时记 nan。"""
    f = e5_file(tag, f"control_n2_T{theta:g}_{case}_tracks.csv")
    if not os.path.exists(f):
        return None
    rows = io_utils.read_csv(f)
    z = np.array([fnum(r["z"]) for r in rows])
    x = np.array([fnum(r["x"]) for r in rows])
    y = np.array([fnum(r["y"]) for r in rows])
    c = np.array([fnum(r["charge"]) for r in rows])
    zs, ang, sep = np.unique(z), [], []
    for zz in zs:
        s = z == zz
        if s.sum() != 2 or (c[s] < 0).any():
            ang.append(np.nan)
            sep.append(np.nan)
            continue
        dx, dy = np.diff(x[s])[0], np.diff(y[s])[0]
        ang.append(np.arctan2(dy, dx))
        sep.append(np.hypot(dx, dy))
    ang, sep = np.array(ang), np.array(sep)
    ok = np.isfinite(ang)
    first_bad = float(zs[~ok][0]) if (~ok).any() else None
    unwrapped = np.full_like(ang, np.nan)
    keep = ok & (zs < (first_bad if first_bad is not None else np.inf))
    unwrapped[keep] = np.unwrap(2 * ang[keep]) / 2          # 轴是无向的: 对 2*phi 展开
    return dict(z=zs, axis=unwrapped, sep=np.where(keep, sep, np.nan), breakdown=first_bad)


def free_rotation_rate(b_disk):
    """圆盘 (R=10) m=2 涡旋在 b_disk 处 L=-2 (Krein +) 最低支频率的一半 (由 rho=8 的通道谱按标度律换算并二次外推)。"""
    f = os.path.join(ROOT, "E3_stability", "hq_hib", "e3b_disk_m2_channels.csv")
    if not os.path.exists(f):
        return float("nan"), float("nan")
    pts = sorted({(round(fnum(r["b"]), 6), fnum(r["lam_im"])) for r in io_utils.read_csv(f)
                  if r["L"] == "-2" and r["krein"] == "1" and fnum(r["lam_im"]) < 0.6})
    pts = np.array(pts)
    s = (R_SAMPLE / RHO_CHANNELS) ** 2
    bp = b_disk * s
    fit = np.polyfit(pts[-5:, 0], pts[-5:, 1], 2)
    w = float(np.polyval(fit, bp)) / s
    return w, 0.5 * w


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    summary = {}
    style.apply()
    fig, axs = style.plt.subplots(2, 3, figsize=(13.5, 7.6))

    # ---- (a) r_s vs Theta
    ax = axs[0, 0]
    rows = [r for r in io_utils.read_csv(e5_file(a.tag, "heating_n2_m2.csv")) if fnum(r["Theta"]) > 0]
    th = np.array([fnum(r["Theta"]) for r in rows])
    r_det = np.array([fnum(r["satellite_r_mean"]) for r in rows])
    r_fit = np.array([fnum(r["r_s_fit"]) for r in rows])
    o = np.argsort(th)
    th, r_det, r_fit = th[o], r_det[o], r_fit[o]
    det = r_det > 0
    ax.loglog(th, r_fit, "o", color=style.PALETTE[0], mfc="none", mec=style.PALETTE[0], ms=6,
              label="harmonic fit of the core")
    ax.loglog(th[det], r_det[det], "s", color=style.PALETTE[1], mec=style.PALETTE[1], ms=5,
              label="detected phase singularities")
    slope = np.polyfit(np.log(th[:4]), np.log(r_fit[:4]), 1)[0]
    tt = np.array([th.min(), th.max()])
    ax.loglog(tt, r_fit[0] * (tt / th[0]) ** 0.5, "-", color=style.INK_SECONDARY, lw=1.0,
              label=f"Theta^(1/2)  (fitted slope at small Theta: {slope:.3f})")
    ax.axhline(2 * H_GRID, color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax.text(th.min() * 1.1, 2 * H_GRID * 1.08, "2h: resolution limit", fontsize=7, color=style.INK_SECONDARY)
    ks = [int(np.argmin(np.abs(th - T))) for T in (0.5, 1.0)]
    ax.loglog(th[ks], r_det[ks], "o", ms=12, mfc="none", mec=style.INK, mew=1.0,
              label="start states of the E5b runs (Theta = 0.5, 1)")
    ax.set_xlabel("quadrupole heating amplitude  Theta")
    ax.set_ylabel("half separation of the twin vortex  r_s")
    ax.set_title("(a) m=2 splits into two unit vortices, r_s ~ Theta^(1/2)", loc="left")
    ax.legend(loc="lower right", fontsize=7)
    summary.update(rs_slope_small_theta=slope, rs_at_theta_0p5=float(r_det[np.argmin(np.abs(th - 0.5))]),
                   rs_at_theta_1=float(r_det[np.argmin(np.abs(th - 1.0))]))

    # ---- (b) (Theta, b) 失稳窗口
    ax = axs[0, 1]
    win = {}
    for r in io_utils.read_csv(e5_file(a.tag, "e5c_windows.csv")):
        if r["n"] == "2" and r["m"] == "2":
            win.setdefault(fnum(r["Theta"]), r)
    summ = {}
    for r in io_utils.read_csv(e5_file(a.tag, "e5c_summary.csv")):
        if r["n"] == "2" and r["m"] == "2":
            summ.setdefault(fnum(r["Theta"]), r)
    Ts = sorted(summ)
    first = {True: True, False: True}
    for T in Ts:
        r = win.get(T)
        if r is not None:
            osc = abs(fnum(r["lam_im_at_max"])) > 1e-6
            st = style.series_style(0 if osc else 1)
            lab = ("oscillatory instability (disk window)" if osc else "real instability (anisotropy-driven)")
            ax.plot([T, T], [fnum(r["b_lo"]), fnum(r["b_hi"])], "-", lw=4.5, color=st["color"],
                    solid_capstyle="butt", label=lab if first[osc] else None)
            first[osc] = False
    bU = np.array([fnum(summ[T]["b_at_U_mark"]) for T in Ts])
    gU = np.array([fnum(summ[T]["growth_at_U_mark"]) for T in Ts])
    ax.plot(Ts, bU, "*", ms=10, color=style.PALETTE[2], mec=style.PALETTE[2],
            label=f"U=15 working point (growth <= {np.nanmax(gU):.0e}: stable)")
    hi = np.array([fnum(win[T]["b_hi"]) if T in win else np.nan for T in Ts])
    big = np.array(Ts) >= 0.5
    fit = np.polyfit(np.array(Ts)[big], hi[big], 1)
    theta_c = (np.mean(bU[-2:]) - fit[1]) / fit[0]
    tt = np.linspace(0.5, theta_c + 0.05, 20)
    ax.plot(tt, np.polyval(fit, tt), "--", color=style.PALETTE[1], lw=1.0)
    ax.annotate(f"upper edge, linear\nextrapolation:\nTheta_c ~ {theta_c:.2f}", (theta_c, np.mean(bU[-2:])),
                (1.04, 0.55), fontsize=7, arrowprops=dict(arrowstyle="->", color=style.INK_SECONDARY, lw=0.7))
    ax.text(0.02, 0.02, "bars start at the lowest b computed (existence limit of the family)",
            transform=ax.transAxes, fontsize=6.5, color=style.INK_SECONDARY)
    ax.set_xlim(-0.05, 1.35)
    ax.set_ylim(-0.12, 3.2)
    ax.set_xlabel("quadrupole heating amplitude  Theta")
    ax.set_ylabel("propagation constant  b")
    ax.set_title("(b) linear stability (bars: unstable range of b)", loc="left")
    ax.legend(loc="upper left", fontsize=7)
    summary.update(theta_c_U15_extrapolated=theta_c, upper_edge_slope=fit[0], max_growth_at_U15=float(np.nanmax(gU)))

    # ---- (c) LG 射入后的形成
    ax = axs[0, 2]
    e4 = os.path.join(ROOT, "E4_dynamic_transmutation", a.tag)
    s4 = {r["run"]: r for r in io_utils.read_csv(os.path.join(e4, "e4b_summary.csv"))} \
        if os.path.exists(os.path.join(e4, "e4b_summary.csv")) else {}
    # 1% 噪声的发射 (e4b_heat2_m2_U15_noise0.01) 形成距离相同 (xi_t50 = 4), 但噪声使谐波拟合的 r_s 抖动, 不画
    for i, (run, lab) in enumerate((("e4b_heat2_m2_U15", "LG m=2 launch"),)):
        f = os.path.join(e4, f"{run}_timeseries.csv")
        if not os.path.exists(f):
            continue
        d = io_utils.read_csv(f)
        z = np.array([fnum(r["z"]) for r in d])
        rs = np.array([fnum(r["r_s_fit"]) for r in d])
        sel = z <= 60
        st = style.series_style(i)
        ax.plot(z[sel], rs[sel], "-", color=st["color"], lw=1.1, label=lab)
    r4 = s4.get("e4b_heat2_m2_U15")
    if r4 is not None:
        ax.axhline(fnum(r4["r_s_stationary"]), color=style.INK_SECONDARY, lw=0.9, ls="--")
        ax.axvline(fnum(r4["xi_t50"]), color=style.INK_SECONDARY, lw=0.8, ls=":")
        ax.text(0.98, 0.97, f"dashed: stationary r_s* = {fnum(r4['r_s_stationary']):.2f}\n"
                f"dotted: xi_t50 = {fnum(r4['xi_t50']):g}\n"
                f"C_2 kept to xi=600 (non-mu weight <= {fnum(r4['nonmu_max']):.0e})",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color=style.INK_SECONDARY)
        summary.update(launch_xi_t50=fnum(r4["xi_t50"]), launch_rs_plateau=fnum(r4["r_s_plateau"]),
                       launch_rs_plateau_std=fnum(r4["r_s_plateau_std"]))
    ax.set_ylim(0, 1.75)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("core radius r_s (harmonic fit)")
    ax.set_title("(c) LG m=2 launch into the heated disk, Theta=0.5", loc="left")

    # ---- (d) 旋转跟随: 轴相对图样的滞后
    ax = axs[1, 0]
    k = 0
    for T in (0.5, 1.0):
        for case in ("rotate_slow", "rotate_fast"):
            tr = pair_track(a.tag, T, case)
            if tr is None:
                continue
            Om = OMEGA[case]
            lag = tr["axis"] - (tr["axis"][0] + Om * tr["z"])
            lag = (lag + np.pi / 2) % np.pi - np.pi / 2
            st = style.series_style(k)
            ax.plot(tr["z"], lag, "-", color=st["color"], lw=1.0,
                    label=f"Theta={T:g}, Omega={Om:g}" + (f" (lock lost at xi={tr['breakdown']:.0f})"
                                                           if tr["breakdown"] is not None else ""))
            if tr["breakdown"] is not None:
                zb = tr["breakdown"]
                kk = int(np.nanargmin(np.abs(np.where(np.isfinite(lag), tr["z"], np.inf) - zb)))
                ax.plot(tr["z"][kk], lag[kk], "X", ms=9, color=st["color"], mec=st["color"])
            summary[f"lag_max_T{T:g}_{case}"] = float(np.nanmax(np.abs(lag)))
            summary[f"breakdown_T{T:g}_{case}"] = tr["breakdown"] if tr["breakdown"] is not None else ""
            k += 1
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("pair axis minus pattern angle (rad)")
    ax.set_title("(d) rotating pattern: the pair follows", loc="left")
    ax.legend(loc="lower left", fontsize=6.8)

    # ---- (e) 关断后的自由共转
    ax = axs[1, 1]
    zs_off = 150.0
    b0 = float("nan")
    disk_row = [r for r in io_utils.read_csv(e5_file(a.tag, "e5c_summary.csv"))
                if r["n"] == "2" and fnum(r["Theta"]) == 0.0]
    if disk_row:
        b0 = fnum(disk_row[0]["b_at_U_mark"])
    w_mode, rate_pred = free_rotation_rate(b0)
    rates = {}
    for i, T in enumerate((0.5, 1.0)):
        tr = pair_track(a.tag, T, "switch_off")
        if tr is None:
            continue
        rot = tr["axis"] - tr["axis"][0]
        st = style.series_style(i)
        ax.plot(tr["z"], rot, "-", color=st["color"], lw=1.3, label=f"Theta={T:g}")
        late = tr["z"] >= zs_off + 50
        rates[T] = float(np.polyfit(tr["z"][late], rot[late], 1)[0])
    ax.axvline(zs_off, color=style.INK_SECONDARY, lw=0.8, ls=":")
    ax.text(zs_off - 6, 25, "heating\nswitched off", fontsize=7, color=style.INK_SECONDARY, ha="right")
    if np.isfinite(rate_pred):
        zz = np.array([zs_off, 600.0])
        ax.plot(zz, rate_pred * (zz - zs_off), "--", color=style.INK, lw=1.0,
                label=f"omega(L=-2)/2 = {rate_pred:.4f} rad per unit xi")
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("rotation of the pair axis (rad)")
    ax.set_title("(e) heating switched off: free co-rotation", loc="left")
    ax.legend(loc="upper left", fontsize=7)
    summary.update(b_disk_U15=b0, omega_L_minus2=w_mode, corotation_predicted=rate_pred,
                   corotation_T0p5=rates.get(0.5, float("nan")), corotation_T1=rates.get(1.0, float("nan")))

    # ---- (f) 间距
    ax = axs[1, 2]
    k = 0
    for case, lab in (("rotate_slow", "slow rotation"), ("rotate_fast", "fast rotation"),
                      ("switch_off", "switched off at xi=150")):
        tr = pair_track(a.tag, 0.5, case)
        if tr is None:
            continue
        st = style.series_style(k)
        sel = (tr["z"] >= 100) & (tr["z"] <= 260)
        ax.plot(tr["z"][sel], tr["sep"][sel], "-", color=st["color"], lw=1.2, label=lab)
        summary[f"sep_range_T0.5_{case}"] = f"{np.nanmin(tr['sep']):.3f}-{np.nanmax(tr['sep']):.3f}"
        k += 1
    ax.axvline(zs_off, color=style.INK_SECONDARY, lw=0.8, ls=":")
    ax.set_ylim(0, 2.6)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("separation of the two vortices  d")
    ax.set_title("(f) separation at Theta=0.5 (zoom on the switch-off)", loc="left")
    ax.legend(loc="lower right", fontsize=7)

    io_utils.write_csv(os.path.join(out, "twin_vortex_summary.csv"), [summary])
    style.save(fig, os.path.join(out, "twin_vortex.png"))
    print("twin vortex summary:")
    for key, v in summary.items():
        print(f"  {key}: {v}")
    print(f"written to {out}")


if __name__ == "__main__":
    main()
