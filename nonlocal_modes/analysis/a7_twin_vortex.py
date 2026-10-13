"""
a7_twin_vortex.py —— 四极 (C_2) 加热下的 "双涡旋": 可调、稳定、可旋转、可释放 (论文 Fig.6 + SI 图)
==============================================================================================
只读取已上传的结果表 (圆盘 R=10, h=0.1, U=15):
  E5a  results/E5_thermal_reconfiguration/<tag>/heating_n2_m2.csv            卫星半径 r_s vs Theta
  E5c  results/E5_thermal_reconfiguration/<tag>/e5c_windows.csv, e5c_summary.csv   (2:2) 的失稳窗口与 U=15 处的增长率
  E5d  results/E5_thermal_reconfiguration/<tag>/e5d_n2_m2_U15_summary.csv     U=15 低频谱: 转动/分裂模 omega_rot、
                                                                              L=-4 模 omega_4、破坏 C_2 的软模 omega_odd
  E5b  results/E5_thermal_reconfiguration/<tag>_omega/control_runs.csv       Omega 扫描 (斜坡启动) 的稳健跟随指标
  E5e  results/E5_thermal_reconfiguration/<tag>/e5e_n2_m2_U15_summary.csv, _T<Theta>_branch.csv  旋转坐标系中的同步定态 (第五轮, 可缺)
       results/E5_thermal_reconfiguration/<tag>/control_n2_T{0.5,1}_*_tracks.csv   第三轮的关断与突然启动的旋转
  E4b  results/E4_dynamic_transmutation/<tag>/e4b_heat2_m2_U15_timeseries.csv  LG m=2 射入加热样品后双涡旋的形成

Fig.6 (twin_vortex.png):
  (a) r_s vs Theta: 渐近 r_s ∝ Theta^(1/2) (牛顿多边形: n=2, m=2)
  (b) (Theta, b) 稳定性: 失稳窗口 + U=15 工作点 (稳定 / 失稳) + E5d 的 Theta_c
  (c) U=15 处破坏 C_2 的平移模: omega_odd^2 (失稳后为 -gamma^2) 对 Theta 线性 -> 叉式分岔 Theta_c
  (d) 跟随相图 (Omega, Theta): 锁定 / 失锁, 与受迫分裂模的一阶共振 omega_rot/2、二阶共振 omega_4/4 对照;
      有 E5e 结果时加上旋转坐标系同步定态的稳定段 / 失稳段 (理论边界, 半透明横条)
  (e) 锁定时卫星半径平方的放大 r2_ratio 对零参数理论 omega_rot / (omega_rot - 2 Omega); 有 E5e 时加同步定态的 r2 (点线)
  (f) 关断后的自由共转: 角速度 = 圆盘 L=-2 分裂模频率的一半
SI (twin_vortex_SI.png): LG 发射中的形成、第三轮突然启动的取向滞后、关断前后的间距
场图 (twin_vortex_fields.png, Fig.6 顶行): E5d --fields_only 导出的 U=15 定态 (results/E5_thermal_reconfiguration/<tag>/fields/),
  每个 Theta 一格: 色相 = 相位、亮度 = sqrt 振幅, 白圈 = +1 奇点, 标注中心荷、两个 +1 涡旋的半间距与稳定性 (E5d)
输出: results/analysis/<tag>/twin_vortex.png、twin_vortex_SI.png、twin_vortex_fields.png (有场文件时)、twin_vortex_summary.csv
运行:  python a7_twin_vortex.py                 (全部)
       python a7_twin_vortex.py --fields_only   (只画场图)
"""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import fieldplot  # noqa: E402
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


def e5e_intervals(tag):
    """E5e 汇总表: {Theta: dict(stable=[(lo, hi, 段)], unstable=[(lo, hi, 段)])}; 段 0 = 从静态双涡旋连续延拓的同相支,
    其余 = 共振另一侧的反相支。没有时返回 {}。"""
    import re
    f = e5_file(tag, "e5e_n2_m2_U15_summary.csv")
    if not os.path.exists(f):
        return {}
    pat = r"([+-]?\d*\.?\d+)\.\.([+-]?\d*\.?\d+)(?:\[([+-]?\d+)\])?"
    out = {}
    for r in io_utils.read_csv(f):
        iv = {k: [(float(a), float(b), int(c or 0)) for a, b, c in re.findall(pat, r.get(f"{k}_intervals", "") or "")]
              for k in ("stable", "unstable")}
        out[fnum(r["Theta"])] = iv
    return out


def e5e_branch(tag, T):
    """E5e 某个 Theta 的同步定态: (Omega, r2_ratio, 是否稳定或未算谱) 数组; 没有时返回 None。"""
    f = e5_file(tag, f"e5e_n2_m2_U15_T{T:g}_branch.csv")
    if not os.path.exists(f):
        return None
    rows = sorted(io_utils.read_csv(f), key=lambda r: fnum(r["Omega"]))
    return (np.array([fnum(r["Omega"]) for r in rows]), np.array([fnum(r["r2_ratio"]) for r in rows]),
            np.array([fnum(r["segment"]) for r in rows]))


def e5d_summary(tag):
    f = e5_file(tag, "e5d_n2_m2_U15_summary.csv")
    if not os.path.exists(f):
        return None
    rows = sorted(io_utils.read_csv(f), key=lambda r: fnum(r["Theta"]))
    return {k: np.array([fnum(r[k]) for r in rows]) for k in ("Theta", "omega_m2", "omega_m4", "omega_odd",
                                                                "max_growth", "Theta_c_fit")}


def twin_fields(tag, out, E, half=3.0):
    """Fig.6 顶行: 各 Theta 的 U=15 定态场图 (E5d --fields_only 的导出)。返回 {Theta: 半间距}。"""
    files = glob.glob(e5_file(tag, os.path.join("fields", "e5d_n2_m2_U15_T*.npz")))
    if not files:
        print(f"no E5d fields in {e5_file(tag, 'fields')}: run e5d_twin_modes.py --fields_only first")
        return {}
    loaded = sorted((fieldplot.load_field(f) for f in files), key=lambda t: t[2]["Theta"])
    k = len(loaded)
    fig = style.plt.figure(figsize=(2.3 * k + 0.3, 3.2))
    gs = fig.add_gridspec(2, k, height_ratios=[1, 0.05], hspace=0.3, wspace=0.08, left=0.02, right=0.98,
                          top=0.8, bottom=0.12)
    half_sep = {}
    for j, (g, q, meta) in enumerate(loaded):
        T = float(meta["Theta"])
        ax = fig.add_subplot(gs[0, j])
        sings, central = fieldplot.draw_field(ax, g, q, half=half)
        core = sorted((np.hypot(x, y), x, y) for (x, y, c) in sings if c > 0 and np.hypot(x, y) < 1.5)
        lines = [f"central charge {central:+d}" if central else "central charge 0"]
        if len(core) == 2:
            rs = 0.5 * np.hypot(core[0][1] - core[1][1], core[0][2] - core[1][2])
            half_sep[T] = rs
            lines.append(f"two +1 vortices, r_s = {rs:.2f}")
        grow = float("nan")
        if E is not None and np.any(np.isclose(E["Theta"], T)):
            grow = float(E["max_growth"][np.isclose(E["Theta"], T)][0])
        ax.set_title(f"Theta = {T:g}", fontsize=9.5)
        if np.isfinite(grow):
            lines.append("unstable, gamma = %.3f" % grow if grow > 1e-6 else "stable")
        ax.text(0.03, 0.03, "\n".join(lines), transform=ax.transAxes, ha="left", va="bottom", fontsize=6.6,
                color="white")
    axb = fig.add_subplot(gs[1, (k - 1) // 2])
    fieldplot.phase_bar(axb)
    fig.suptitle(f"quadrupole-heated disk (R = 10), U = 15: hue = phase, brightness = sqrt(amplitude)\n"
                 f"white circles: +1 vortices; window |x|, |y| <= {half:g}", fontsize=8.5)
    style.save(fig, os.path.join(out, "twin_vortex_fields.png"), tight=False)
    return half_sep


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    p.add_argument("--fields_only", action="store_true", help="只画场图 (twin_vortex_fields.png)")
    p.add_argument("--e5e_tag", default=None, help="E5e 结果所在的 tag (默认与 --tag 相同; 云端检查用 cloudtest_e5e)")
    a = p.parse_args()
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    summary = {}
    style.apply()
    if a.fields_only:
        for T, rs in twin_fields(a.tag, out, e5d_summary(a.tag)).items():
            print(f"  Theta={T:g}: r_s from the field = {rs:.4f}")
        print(f"written to {out}")
        return
    fig, axs = style.plt.subplots(2, 3, figsize=(14.0, 8.0))
    E = e5d_summary(a.tag)

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
    ax.set_xlabel("quadrupole heating amplitude  Theta")
    ax.set_ylabel("half separation of the twin vortex  r_s")
    ax.set_title("(a) m=2 splits into two unit vortices", loc="left")
    ax.legend(loc="lower right", fontsize=7)
    summary.update(rs_slope_small_theta=slope)

    # ---- (b) (Theta, b) 稳定性
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
            # Im λ = 0 的窗口: 低 b 端是 κ=0 (L=-2) 模, Theta >= 0.75 时上端由破坏 C_2 的 κ=1 平移模决定 (Progress.md C10)
            lab = ("oscillatory instability (disk window)" if osc else "non-oscillatory instability (Im lambda = 0)")
            ax.plot([T, T], [fnum(r["b_lo"]), fnum(r["b_hi"])], "-", lw=4.5, color=st["color"],
                    solid_capstyle="butt", label=lab if first[osc] else None)
            first[osc] = False
    bU = np.array([fnum(summ[T]["b_at_U_mark"]) for T in Ts])
    unst = np.array([summ[T].get("unstable_at_U_mark", "") == "True" for T in Ts])
    Ta = np.array(Ts)
    ax.plot(Ta[~unst], bU[~unst], "*", ms=10, color=style.PALETTE[2], mec=style.PALETTE[2],
            label="U=15 working point: stable")
    if unst.any():
        ax.plot(Ta[unst], bU[unst], "X", ms=9, color=style.PALETTE[7], mec=style.PALETTE[7],
                label="U=15 working point: unstable (kappa=1)")
    if E is not None and np.isfinite(E["Theta_c_fit"][0]):
        tc = float(E["Theta_c_fit"][0])
        ax.axvline(tc, color=style.INK_SECONDARY, lw=0.9, ls="--")
        ax.text(tc + 0.02, 0.15, f"Theta_c = {tc:.2f}\n(panel c)", fontsize=7, color=style.INK_SECONDARY)
        summary["theta_c_e5d"] = tc
    ax.text(0.02, 0.02, "bars start at the lowest b computed (existence limit)", transform=ax.transAxes,
            fontsize=6.5, color=style.INK_SECONDARY)
    ax.set_xlim(-0.05, 1.6)
    ax.set_ylim(-0.12, 3.25)
    ax.set_xlabel("quadrupole heating amplitude  Theta")
    ax.set_ylabel("propagation constant  b")
    ax.set_title("(b) linear stability (bars: unstable b)", loc="left")
    ax.legend(loc="upper left", fontsize=6.8)

    # ---- (c) C_2 破缺的平移模软化 -> 叉式分岔
    ax = axs[0, 2]
    if E is not None:
        st_ok = (E["max_growth"] < 1e-5) & np.isfinite(E["omega_odd"]) & (E["omega_odd"] > 0)
        w2 = np.where(st_ok, E["omega_odd"] ** 2, -E["max_growth"] ** 2)
        ax.plot(E["Theta"][st_ok], w2[st_ok], "o", ms=6, color=style.PALETTE[0], mec=style.PALETTE[0],
                label="stable: omega_odd^2 (E5d)")
        ax.plot(E["Theta"][~st_ok], w2[~st_ok], "s", ms=6, color=style.PALETTE[7], mec=style.PALETTE[7],
                label="unstable: -gamma^2 (E5d)")
        gE = np.array([fnum(summ[T]["growth_at_U_mark"]) if summ[T].get("unstable_at_U_mark") == "True" else np.nan
                       for T in Ts])
        okE = np.isfinite(gE)
        ax.plot(Ta[okE], -gE[okE] ** 2, "D", ms=9, mfc="none", color=style.PALETTE[1], mec=style.PALETTE[1],
                label="unstable: -gamma^2 (E5c, independent run)")
        fit = np.polyfit(E["Theta"][st_ok], w2[st_ok], 1)
        tt = np.linspace(0, E["Theta"].max() + 0.05, 50)
        ax.plot(tt, np.polyval(fit, tt), "-", color=style.INK_SECONDARY, lw=1.0,
                label=f"linear fit of the stable points: {fit[1]:.4f} - {-fit[0]:.4f} Theta")
        ax.axhline(0, color=style.INK_SECONDARY, lw=0.7)
        summary.update(odd_w2_intercept=fit[1], odd_w2_slope=fit[0])
    ax.set_xlabel("quadrupole heating amplitude  Theta")
    ax.set_ylabel("omega_odd^2  (or -gamma^2)   at U=15")
    ax.set_title("(c) C_2-breaking translation mode softens", loc="left")
    ax.legend(loc="lower left", fontsize=6.8)

    # ---- (d) 跟随相图
    ax = axs[1, 0]
    scan = e5_file(f"{a.tag}_omega", "control_runs.csv")
    runs = [r for r in io_utils.read_csv(scan) if r["case"].startswith("rotate_Om")] if os.path.exists(scan) else []
    for lab, sel, mk, col in (("locked (rho = 1, core intact)", True, "o", style.PALETTE[0]),
                              ("lock lost (phase slip or core change)", False, "X", style.PALETTE[1])):
        pts = [(fnum(r["omega"]), fnum(r["Theta"])) for r in runs if (r["locked"] == "True") == sel]
        if pts:
            ax.plot([q[0] for q in pts], [q[1] for q in pts], mk, ms=7 if mk == "o" else 8, color=col, mec=col,
                    ls="none", label=lab)
    # E5e: 旋转坐标系中同步定态的稳定段 / 失稳段 (理论边界; 第五轮)
    EE = e5e_intervals(a.e5e_tag or a.tag)
    styles = {("stable", True): (style.PALETTE[0], 7, 0.30, "co-rotating state stable (in-phase branch, E5e)"),
              ("stable", False): (style.PALETTE[2], 7, 0.30, "stable beyond a resonance (anti-phase branch, E5e)"),
              ("unstable", True): (style.PALETTE[1], 3, 0.55, "co-rotating state unstable (E5e)"),
              ("unstable", False): (style.PALETTE[1], 3, 0.55, None)}
    used = set()
    for T, iv in sorted(EE.items()):
        for kind in ("unstable", "stable"):
            for lo, hi, sg in iv[kind]:
                col, lw, al, lab = styles[(kind, sg == 0)]
                dy = 0.0 if sg == 0 else (0.035 if kind == "stable" else 0.035)
                ax.plot([lo, hi], [T + dy, T + dy], "-", color=col, lw=lw, alpha=al, solid_capstyle="butt", zorder=1,
                        label=None if (lab is None or lab in used) else lab)
                if lab:
                    used.add(lab)
    if E is not None:
        ok = E["max_growth"] < 1e-5
        ax.plot(E["omega_m2"][ok] / 2, E["Theta"][ok], "-", color=style.INK, lw=1.1,
                label="1st-order resonance omega_rot/2")
        ax.plot(E["omega_m4"][ok] / 4, E["Theta"][ok], "--", color=style.INK, lw=1.1,
                label="2nd-order resonance omega_4/4 (coupling ~ Theta^2)")
        # 反向旋转的失锁 (Theta >= 0.75) 与 C_2 破缺的软模 (omega_odd <= 0.09) 同量级, 但 Theta <= 0.5 时 |Omega| > omega_odd
        # 仍锁定 -> 没有简单的 |Omega| = omega_odd 判据, 不画线 (Progress.md C10)
    ax.set_xlim(-0.22, 0.17)
    ax.set_ylim(0, 1.62)
    ax.set_xlabel("pattern angular velocity  Omega   (>0: co-rotating with the free pair)")
    ax.set_ylabel("quadrupole heating amplitude  Theta")
    ax.set_title("(d) following map of the rotating pattern (ramped)", loc="left")
    ax.legend(loc="upper left", fontsize=6.8, ncol=1)
    if runs:
        for T in sorted({fnum(r["Theta"]) for r in runs}):
            lk = [fnum(r["omega"]) for r in runs if fnum(r["Theta"]) == T and r["locked"] == "True"]
            summary[f"locked_Omega_range_T{T:g}"] = f"{min(lk):g}..{max(lk):g}" if lk else ""

    # ---- (e) r2_ratio 对零参数理论
    ax = axs[1, 1]
    for i, T in enumerate(sorted({fnum(r["Theta"]) for r in runs})):
        st = style.series_style(i)
        mine = [r for r in runs if fnum(r["Theta"]) == T and r["r2_ratio"] not in ("", None)]
        lk = [r for r in mine if r["locked"] == "True"]
        nl = [r for r in mine if r["locked"] != "True"]
        ax.plot([fnum(r["omega"]) for r in lk], [fnum(r["r2_ratio"]) for r in lk], st["marker"], ms=6,
                color=st["color"], mec=st["color"], ls="none", label=f"Theta={T:g} (locked)")
        if nl:
            ax.plot([fnum(r["omega"]) for r in nl], [fnum(r["r2_ratio"]) for r in nl], st["marker"], ms=6,
                    color=st["color"], mfc="none", mec=st["color"], ls="none")
        br = e5e_branch(a.e5e_tag or a.tag, T)
        if br is not None:
            # 只画稳定的同步定态 (E5e 汇总表中该段的稳定区间内)
            stab = EE.get(T, {}).get("stable", [])
            for sg in np.unique(br[2]):
                ok = np.zeros(len(br[0]), bool)
                for lo, hi, s2 in stab:
                    if s2 == sg:
                        ok |= (br[0] >= lo - 1e-9) & (br[0] <= hi + 1e-9)
                sel = (br[2] == sg) & np.isfinite(br[1]) & ok
                ax.plot(br[0][sel], br[1][sel], ":", color=st["color"], lw=1.4)
        if E is not None:
            w = float(np.interp(T, E["Theta"], E["omega_m2"]))
            oo = np.linspace(-0.21, 0.92 * w / 2, 120)
            ax.plot(oo, w / (w - 2 * oo), "-", color=st["color"], lw=1.0)
            dev = [fnum(r["r2_ratio"]) / (w / (w - 2 * fnum(r["omega"]))) - 1 for r in lk
                   if abs(fnum(r["omega"])) <= 0.1 and abs(fnum(r["omega"]) - 0.025) > 0.01]
            if dev:
                summary[f"r2_theory_maxdev_T{T:g}"] = float(np.max(np.abs(dev)))
    ax.set_yscale("log")
    ax.set_xlabel("pattern angular velocity  Omega")
    ax.set_ylabel("<r_s^2> / r_s0^2   after the ramp")
    ax.set_title("(e) separation response: omega_rot / (omega_rot - 2 Omega)", loc="left")
    ax.text(0.98, 0.03, "solid lines: no free parameter (omega_rot from E5d)\ndotted: stable co-rotating stationary state (E5e)"
            "\nopen symbols: lock lost", transform=ax.transAxes, ha="right", fontsize=6.8, color=style.INK_SECONDARY)
    ax.legend(loc="upper left", fontsize=6.8)

    # ---- (f) 关断后的自由共转
    ax = axs[1, 2]
    zs_off = 150.0
    if E is not None and E["Theta"][0] == 0:
        w_mode = float(E["omega_m2"][0])
    else:
        disk_row = [r for r in io_utils.read_csv(e5_file(a.tag, "e5c_summary.csv"))
                    if r["n"] == "2" and fnum(r["Theta"]) == 0.0]
        w_mode = free_rotation_rate(fnum(disk_row[0]["b_at_U_mark"]))[0] if disk_row else float("nan")
    rate_pred = 0.5 * w_mode
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
    ax.set_title("(f) heating switched off: free co-rotation", loc="left")
    ax.legend(loc="upper left", fontsize=7)
    summary.update(omega_L_minus2=w_mode, corotation_predicted=rate_pred,
                   corotation_T0p5=rates.get(0.5, float("nan")), corotation_T1=rates.get(1.0, float("nan")))
    style.save(fig, os.path.join(out, "twin_vortex.png"))

    # ================= SI: 形成、第三轮突然启动的滞后、间距
    fig, axs = style.plt.subplots(1, 3, figsize=(13.5, 3.8))
    ax = axs[0]
    e4 = os.path.join(ROOT, "E4_dynamic_transmutation", a.tag)
    s4 = {r["run"]: r for r in io_utils.read_csv(os.path.join(e4, "e4b_summary.csv"))} \
        if os.path.exists(os.path.join(e4, "e4b_summary.csv")) else {}
    f4 = os.path.join(e4, "e4b_heat2_m2_U15_timeseries.csv")
    if os.path.exists(f4):
        d = io_utils.read_csv(f4)
        z = np.array([fnum(r["z"]) for r in d])
        rs = np.array([fnum(r["r_s_fit"]) for r in d])
        sel = z <= 60
        ax.plot(z[sel], rs[sel], "-", color=style.PALETTE[0], lw=1.1)
    r4 = s4.get("e4b_heat2_m2_U15")
    if r4 is not None:
        ax.axhline(fnum(r4["r_s_stationary"]), color=style.INK_SECONDARY, lw=0.9, ls="--")
        ax.axvline(fnum(r4["xi_t50"]), color=style.INK_SECONDARY, lw=0.8, ls=":")
        ax.text(0.98, 0.97, f"dashed: stationary r_s* = {fnum(r4['r_s_stationary']):.2f}\n"
                f"dotted: xi_t50 = {fnum(r4['xi_t50']):g}\n"
                f"C_2 kept to xi=600 (non-mu weight <= {fnum(r4['nonmu_max']):.0e})",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color=style.INK_SECONDARY)
        summary.update(launch_xi_t50=fnum(r4["xi_t50"]))
    ax.set_ylim(0, 1.75)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("core radius r_s (harmonic fit)")
    ax.set_title("(a) LG m=2 launch into the heated disk, Theta=0.5", loc="left")

    ax = axs[1]
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
                    label=f"Theta={T:g}, Omega={Om:g}" + (" (core change)" if tr["breakdown"] is not None else ""))
            k += 1
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("pair axis minus pattern angle (rad)")
    ax.set_title("(b) sudden start of the rotation (round 3)", loc="left")
    ax.legend(loc="lower left", fontsize=6.8)

    ax = axs[2]
    k = 0
    for case, lab in (("rotate_slow", "slow rotation"), ("rotate_fast", "fast rotation"),
                      ("switch_off", "switched off at xi=150")):
        tr = pair_track(a.tag, 0.5, case)
        if tr is None:
            continue
        st = style.series_style(k)
        sel = (tr["z"] >= 100) & (tr["z"] <= 260)
        ax.plot(tr["z"][sel], tr["sep"][sel], "-", color=st["color"], lw=1.2, label=lab)
        k += 1
    ax.axvline(zs_off, color=style.INK_SECONDARY, lw=0.8, ls=":")
    ax.set_ylim(0, 2.6)
    ax.set_xlabel("propagation distance xi")
    ax.set_ylabel("separation of the two vortices  d")
    ax.set_title("(c) separation at Theta=0.5", loc="left")
    ax.legend(loc="lower right", fontsize=7)
    style.save(fig, os.path.join(out, "twin_vortex_SI.png"))

    for T, rs in twin_fields(a.tag, out, E).items():
        summary[f"rs_field_T{T:g}"] = rs
    io_utils.write_csv(os.path.join(out, "twin_vortex_summary.csv"), [summary])
    print("twin vortex summary:")
    for key, v in summary.items():
        print(f"  {key}: {v}")
    print(f"written to {out}")


if __name__ == "__main__":
    main()
