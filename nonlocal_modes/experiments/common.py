"""
common.py —— 各实验脚本共用的工具
=================================
  * build_sample      : 用名字构造样品 ('rect', 'poly3'...'poly6', 'disk', 'surface')
  * robust_family     : 稳健地求一个解族 (先在若干 b 上找到第一个收敛解, 再向上/向下双向延拓)
  * param_continuation: 沿“样品参数”(尺寸、边界温度等) 延拓, 失败时自动二分加密步长
  * Tee               : 把屏幕输出同时写进 log.txt (日志需要上传)
  * standard_args     : 所有脚本共有的命令行参数 (--tag --quick --h ...)
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


def standard_args(description, default_h=0.2, quick_h=0.4):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--tag", default=None,
                   help="结果子目录名 (默认: 正式运行 'main', --quick 时 'quick')")
    p.add_argument("--quick", action="store_true",
                   help="快速冒烟测试: 粗网格 + 少量参数点, 几分钟内跑完, 只用于检查程序能否运行")
    p.add_argument("--h", type=float, default=None, help=f"网格间距 (默认 {default_h}, quick 时 {quick_h})")
    p._nl_defaults = (default_h, quick_h)
    return p


def finalize_args(p, exp_name):
    args = p.parse_args()
    default_h, quick_h = p._nl_defaults
    if args.h is None:
        args.h = quick_h if args.quick else default_h
    if args.tag is None:
        args.tag = "quick" if args.quick else "main"
    args.exp = exp_name
    args.res_dir = io_utils.results_dir(exp_name, args.tag)
    args.data_dir = io_utils.data_dir(exp_name, args.tag)
    sys.stdout = Tee(os.path.join(args.res_dir, "log.txt"))
    print(f"\n===== {exp_name}  tag={args.tag}  h={args.h}  quick={args.quick}  "
          f"{time.strftime('%Y-%m-%d %H:%M:%S')} =====")
    print(f"results -> {args.res_dir}\nraw data -> {args.data_dir}")
    io_utils.write_run_info(args.res_dir, args)
    return args


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
                  callback=None, verbose=True):
    """
    求一个解族: 先在 b_start 附近找初解, 再分别向 b_high 与 b_low 延拓。
    返回 (按 b 升序的解列表, info dict)。
    callback(sol) 对每个收敛解调用一次。
    """
    if b_start is None:
        b_start = 0.5 * (b_low + b_high)
    cands = [b_start]
    for k in (1, 2, 3):
        for s in (+1, -1):
            b = b_start + s * 0.75 * k
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
                       max_depth=5, verbose=True):
    """
    沿样品参数延拓 (例如矩形宽度 d_eta、多边形保角半径 rho、边界温度幅度 Theta)。
    build_model(v) -> NonlocalModel;  solve_fn(model, w_guess, prev_sol) -> sol (dict, 含 converged)
    targets: 按远离 start_value 的顺序排列的目标参数值。
    某一步牛顿失败时, 在上一个成功值与目标值之间插入中点重试, 最多二分 max_depth 层。
    返回 [(v, model, sol), ...] (只含 targets 中成功的点)。
    """
    out = []
    cur_v, cur_M, cur_s = start_value, start_model, start_sol
    for v in targets:
        pending = [v]
        depth = 0
        while pending:
            vt = pending[-1]
            Mt = build_model(vt)
            sol = solve_fn(Mt, transfer_solution(cur_M, cur_s["w"], Mt), cur_s)
            if sol["converged"]:
                cur_v, cur_M, cur_s = vt, Mt, sol
                pending.pop()
                depth = max(depth - 1, 0)
                if vt == v:
                    out.append((v, Mt, sol))
            else:
                depth += 1
                if depth > max_depth:
                    if verbose:
                        print(f"   param continuation stopped between {cur_v:g} and {vt:g}: {sol['message']}")
                    return out
                pending.append(0.5 * (cur_v + vt))
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
