"""
common.py —— 各实验脚本共用的工具
=================================
  * build_sample      : 用名字构造样品 ('rect', 'poly3'...'poly6', 'disk', 'surface')
  * robust_family     : 稳健地求一个解族 (先在若干 b 上找到第一个收敛解, 再向上/向下双向延拓)
  * param_continuation: 沿“样品参数”(尺寸、边界温度等) 延拓, 失败时自动二分加密步长
  * Tee               : 把屏幕输出同时写进 log.txt (日志需要上传)
  * standard_args     : 所有脚本共有的命令行参数 (--tag --quick --hq --h ...)
  * PROFILES / pick   : 精度档位 quick / main / hq 及各档参数
"""

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, io_utils, newton, plotting  # noqa: E402
from nlmodes.geometry import (disk_sample, polygon_sample, rectangle_sample,  # noqa: E402
                              surface_interface_sample)
from nlmodes.model import NonlocalModel  # noqa: E402


class Tee:
    """同时输出到屏幕和文件。"""

    def __init__(self, path):
        self.file = open(path, "a", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, s):
        self.stdout.write(s)
        self.file.write(s)
        self.file.flush()

    def flush(self):
        self.stdout.flush()
        self.file.flush()


# ---------------------------------------------------------------------------
# 精度档位 (precision profiles)
#   quick : 冒烟测试, 粗网格, 几分钟内跑完, 结果不上传
#   main  : 普通正式运行 (h = 0.2)
#   hq    : 高精度正式运行 (--hq), 用于论文数据; 计算量约为 main 的 10-20 倍, 建议在超算上拆分并行
# 说明: 牛顿迭代的收敛判据不能设得比 ~1e-13 更小 —— 双精度浮点的机器精度约 2.2e-16, 乘上差分算子
# 的系数后残差的舍入噪声在 1e-14~1e-13 量级; 设成 1e-20 会使所有求解都“不收敛”。
# 论文结果的误差主要来自网格 h 和步长 dz (1e-3~1e-5 量级), 所以 hq 档把算力花在加密网格与步长上。
# ---------------------------------------------------------------------------
PROFILES = {
    #        光学/热网格 h   大样品网格 h_large  牛顿判据   传播步长   b 步长缩放
    "quick": dict(h=None, h_large=None, tol=1e-9, dz=0.02, db_factor=1.0),
    "main": dict(h=None, h_large=None, tol=1e-9, dz=0.02, db_factor=1.0),
    "hq": dict(h=0.1, h_large=0.15, tol=1e-12, dz=0.01, db_factor=0.5),
}
LARGE_SAMPLES = ("rect", "surface")     # 截面大 (20x40, 40x40) 的样品在 hq 档用 h_large 控制计算量


def standard_args(description, default_h=0.2, quick_h=0.4):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--tag", default=None,
                   help="结果子目录名 (默认: 正式运行 'main', --hq 时 'hq', --quick 时 'quick')")
    p.add_argument("--quick", action="store_true",
                   help="快速冒烟测试: 粗网格 + 少量参数点, 几分钟内跑完, 只用于检查程序能否运行")
    p.add_argument("--hq", action="store_true",
                   help="高精度正式运行: h=0.1 (矩形/界面样品 0.15), 牛顿判据 1e-12, dz=0.01, b 步长减半, "
                        "更长传播与更密的稳定性扫描 (论文数据用这一档)")
    p.add_argument("--job", default=None,
                   help="拆分并行运行时的任务名 (由 hpc/make_jobs.py 自动给出): 日志写入 log_<job>.txt、"
                        "run_info_<job>.json, 避免多个进程同时写同一个文件")
    p.add_argument("--h", type=float, default=None,
                   help=f"网格间距 (默认 {default_h}; quick {quick_h}; hq 0.1, 矩形/界面样品 0.15)。"
                        "显式给出时所有样品都用这个值")
    p.add_argument("--h_scale", type=float, default=1.0,
                   help="把本档的网格间距 (含大样品) 统一乘以此因子, 用于网格收敛性检验, "
                        "例如 --hq --h_scale 1.5 --tag hq_coarse, 再用 compare_tags.py 与 hq 比较")
    p._nl_defaults = (default_h, quick_h)
    return p


def finalize_args(p, exp_name):
    args = p.parse_args()
    if args.quick and args.hq:
        p.error("--quick 与 --hq 不能同时使用")
    default_h, quick_h = p._nl_defaults
    args.profile = "quick" if args.quick else ("hq" if args.hq else "main")
    prof = PROFILES[args.profile]
    if args.h is None:
        args.h = {"quick": quick_h, "main": default_h, "hq": prof["h"]}[args.profile]
        args.h_large = {"quick": quick_h, "main": default_h, "hq": prof["h_large"]}[args.profile]
    else:
        args.h_large = args.h
    if args.h_scale != 1.0:
        args.h = round(args.h * args.h_scale, 6)
        args.h_large = round(args.h_large * args.h_scale, 6)
    args.tol, args.dz, args.db_factor = prof["tol"], prof["dz"], prof["db_factor"]
    newton.DEFAULT_TOL = args.tol
    if args.tag is None:
        args.tag = args.profile
    args.exp = exp_name
    args.res_dir = io_utils.results_dir(exp_name, args.tag)
    args.data_dir = io_utils.data_dir(exp_name, args.tag)
    suffix = f"_{args.job}" if args.job else ""
    sys.stdout = Tee(os.path.join(args.res_dir, f"log{suffix}.txt"))
    print(f"\n===== {exp_name}  tag={args.tag}  profile={args.profile}  h={args.h} (large samples "
          f"{args.h_large})  newton tol={args.tol:g}  dz={args.dz}  {time.strftime('%Y-%m-%d %H:%M:%S')} =====")
    print(f"results -> {args.res_dir}\nraw data -> {args.data_dir}")
    io_utils.write_run_info(args.res_dir, args, filename=f"run_info{suffix}.json")
    return args


def pick(args, quick, main, hq):
    """按精度档位选参数值: pick(args, quick=..., main=..., hq=...)。"""
    return {"quick": quick, "main": main, "hq": hq}[args.profile]


def sample_h(args, kind):
    """该样品使用的网格间距 (hq 档中矩形/界面样品用 h_large)。"""
    return args.h_large if kind in LARGE_SAMPLES else args.h


def build_sample(kind, h, size):
    """
    kind : 'rect' (size=(d_eta, d_zeta)), 'poly<n>' (size=rho, 保角半径), 'disk' (size=R),
           'surface' (size=L)
    返回 (sample, nsym): nsym 为样品的旋转对称阶数 (圆盘为 None, 矩形为 2)。
    """
    if kind == "rect":
        Lx, Ly = size
        return rectangle_sample(Lx, Ly, h), 2
    if kind.startswith("poly"):
        n = int(kind[4:])
        return polygon_sample(n, h, rho=size), n
    if kind == "disk":
        return disk_sample(size, h), None
    if kind == "surface":
        return surface_interface_sample(L=size, h=h), None
    raise ValueError(kind)


def initial_solution(model, make_shape, b_candidates, widths=np.linspace(0.8, 6.0, 14), real=None,
                     verbose=True):
    """在 b_candidates 中依次尝试, 返回第一个收敛的解 (或 None)。"""
    for b in b_candidates:
        try:
            w0, wid = guesses.best_guess(model, make_shape, b, widths)
        except ValueError:
            continue
        sol = newton.solve_robust(model, w0, b, real=real, maxit=30)
        if verbose:
            print(f"  initial solve at b={b:.3f} (guess width {wid:.2f}): {sol['message']}, "
                  f"its={sol['iterations']}")
        if sol["converged"]:
            sol["param"] = b
            return sol
    return None


def robust_family(model, make_shape, b_low, b_high, db=0.25, b_start=None, real=None,
                  callback=None, verbose=True, align_db=False):
    """
    求一个解族: 先在 b_start 附近找初解, 再分别向 b_high 与 b_low 延拓。
    返回 (按 b 升序的解列表, info dict)。
    callback(sol) 对每个收敛解调用一次。
    align_db=True: 起点与备选起点都取在 db 的整数倍上, 于是解族的 b 网格与网格间距 h、与哪个备选起点收敛无关
                   (不同网格的结果可以按 b 逐点比较, 见 compare_tags.py; E2a/E3 使用)。
    """
    if b_start is None:
        b_start = 0.5 * (b_low + b_high)
    off = 0.75
    if align_db:
        step = abs(db)
        b_start = round(b_start / step) * step
        off = max(1, round(0.75 / step)) * step
    cands = [b_start]
    for k in (1, 2, 3):
        for s in (+1, -1):
            b = b_start + s * off * k
            if b_low <= b <= b_high:
                cands.append(b)
    first = initial_solution(model, make_shape, cands, real=real, verbose=verbose)
    if first is None:
        return [], dict(stop_up="no initial solution", stop_down="no initial solution")
    b0 = first["b"]
    up, info_up = newton.continuation(model, first["w"], b0, b_high, abs(db), real=real,
                                      verbose=verbose)
    down, info_down = newton.continuation(model, first["w"], b0, b_low, -abs(db), real=real,
                                          verbose=verbose)
    sols = sorted(down[1:] + up, key=lambda s: s["b"])
    if callback:
        for s in sols:
            callback(s)
    return sols, dict(stop_up=info_up["stop_reason"], stop_down=info_down["stop_reason"],
                      b_min=sols[0]["b"] if sols else None, b_max=sols[-1]["b"] if sols else None)


def pick_closest(sols, key, target):
    """在解列表里挑 key (='b' 或 'power') 最接近 target 的解。"""
    return min(sols, key=lambda s: abs(s[key] - target)) if sols else None


def save_family(model, sols, folder_res, folder_data, name, nsym=None, m=None, plot_bs=None,
                center=(0.0, 0.0), extra_cols=None):
    """
    保存一个解族:
      results/<...>/<name>_summary.csv     每个解一行诊断量 (上传)
      results/<...>/<name>_U_vs_b.png      U(b) 等曲线 (上传)
      results/<...>/<name>_b<val>.png      选定 b 处的模式图 (上传)
      data/<...>/<name>_b<val>.npz         所有解的双精度数据 (不上传, 供后续实验读取)
    """
    rows = []
    for s in sols:
        row = dg.mode_summary(model, s, nsym=nsym, center=center, m=m)
        if extra_cols:
            row.update(extra_cols)
        rows.append(row)
        io_utils.save_solution(os.path.join(folder_data, f"{name}_b{s['b']:.4f}.npz"), model, s,
                               extra=dict(m=m, nsym=nsym))
    io_utils.write_csv(os.path.join(folder_res, f"{name}_summary.csv"), rows)
    if len(rows) >= 2:
        plotting.plot_family(rows, os.path.join(folder_res, f"{name}_U_vs_b.png"), x="b",
                             ys=("U", "n_max", "ring_radius"), title=name)
    for bt in (plot_bs or []):
        s = pick_closest(sols, "b", bt)
        if s is not None and abs(s["b"] - bt) < 0.13:
            plotting.plot_mode(model, s, os.path.join(folder_res, f"{name}_b{s['b']:.2f}.png"),
                               title=f"{name}  b={s['b']:.3f}  U={s['power']:.3f}", center=center)
    return rows


def transfer_solution(model_from, w, model_to):
    """把 model_from 上的解 w 插值到 model_to 的光学点上 (网格相同则直接取值)。"""
    g = model_from.grid
    saved = dict(x=g.x, y=g.y, h=g.h, q=dg.field_grid(model_from, w))
    return io_utils.solution_on_model(saved, model_to)


def param_continuation(build_model, start_value, start_model, start_sol, targets, solve_fn,
                       max_depth=5, verbose=True, secant=True):
    """
    沿样品参数延拓 (例如矩形宽度 d_eta、多边形保角半径 rho、边界温度幅度 Theta)。
    build_model(v) -> NonlocalModel;  solve_fn(model, w_guess, prev_sol) -> sol (dict, 含 converged)
    targets: 按远离 start_value 的顺序排列的目标参数值。
    步长控制: 先一步走到下一个目标; 牛顿失败则步长减半, 之后以减半后的步长逐步逼近该目标 (不再回头重试整步,
    这样在解族端点/折返点附近最多失败约 max_depth 次就停下); 步长小于初始步长 / 2^max_depth 时停止。
    secant=True (默认): 用最近两个解做线性外推作为初值 (割线预测)。解的结构随参数明显变化时
                 (例如加热增强、卫星环外移) 比直接用上一个解作初值好得多, 细网格上尤其如此。
    返回 [(v, model, sol), ...] (只含 targets 中成功的点)。
    """
    out = []
    cur_v, cur_M, cur_s = start_value, start_model, start_sol
    prev = None                                   # 上一个成功点 (v, model, sol), 供割线预测
    for v in targets:
        step = v - cur_v
        min_step = abs(step) / 2 ** max_depth
        while abs(v - cur_v) > 1e-12 * max(1.0, abs(v)):
            vt = v if abs(step) >= abs(v - cur_v) - 1e-12 else cur_v + step
            Mt = build_model(vt)
            w_guess = transfer_solution(cur_M, cur_s["w"], Mt)
            if secant and prev is not None and prev[0] != cur_v:
                w_prev = newton.align_phase(transfer_solution(prev[1], prev[2]["w"], Mt), w_guess)
                w_guess = w_guess + (vt - cur_v) / (cur_v - prev[0]) * (w_guess - w_prev)
            sol = solve_fn(Mt, w_guess, cur_s)
            if sol["converged"]:
                prev = (cur_v, cur_M, cur_s)
                cur_v, cur_M, cur_s = vt, Mt, sol
            else:
                step = step / 2
                if abs(step) < min_step:
                    if verbose:
                        print(f"   param continuation stopped between {cur_v:g} and {vt:g}: {sol['message']}",
                              flush=True)
                    return out
        out.append((v, cur_M, cur_s))
        if verbose:
            print(f"   param continuation: value {v:g} reached (b={cur_s['b']:.4f}, "
                  f"{cur_s.get('iterations', -1)} Newton its)", flush=True)
    return out


def rect_state_via_square(h, m, b, d_eta=20.0, d_zeta=40.0, verbose=True):
    """
    矩形样品中强各向异性的高阶态 (例如 20x40 中的 m>=4) 直接用 LG 初值往往不收敛。
    办法: 先在正方形 d_zeta x d_zeta 中求解 (各向异性弱, LG 初值有效), 再固定 b 沿 d_eta 收缩到目标值
    (参数延拓, 失败时自动二分步长)。返回 (model, sol) 或 (None, None)。
    """
    from nlmodes.geometry import rectangle_sample
    Msq = NonlocalModel(rectangle_sample(d_zeta, d_zeta, h))
    s = initial_solution(Msq, lambda W: guesses.lg_vortex(Msq, m, W), [b, b + 0.75, b - 0.75], verbose=False)
    if s is None:
        return None, None
    if abs(s["b"] - b) > 1e-9:
        sl, _ = newton.continuation(Msq, s["w"], s["b"], b, 0.25 if b > s["b"] else -0.25, verbose=False)
        if not sl or abs(sl[-1]["b"] - b) > 1e-9:
            return None, None
        s = sl[-1]
    steps = list(np.linspace(d_zeta, d_eta, int(round((d_zeta - d_eta) / 4)) + 1)[1:])
    scan = param_continuation(lambda d: NonlocalModel(rectangle_sample(d, d_zeta, h)), d_zeta, Msq, s, steps,
                              lambda model, w0, prev: newton.solve_robust(model, w0, b), verbose=verbose)
    if scan and abs(scan[-1][0] - d_eta) < 1e-9:
        if verbose:
            print(f"   obtained m={m} state in {d_eta:g}x{d_zeta:g} by aspect-ratio continuation from the square")
        return scan[-1][1], scan[-1][2]
    return None, None
