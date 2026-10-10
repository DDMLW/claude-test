"""
io_utils.py —— 数据保存/读取与运行信息记录
==========================================

目录约定 (见 nonlocal_modes/results/README.md):
    nonlocal_modes/data/<实验>/<tag>/...     原始数据 (大, 不上传 GitHub, 已在 .gitignore 中)
    nonlocal_modes/results/<实验>/<tag>/...  整理后的结果 (图、CSV、少量代表性解) —— 需要上传
"""

import csv
import datetime
import json
import os
import platform
import subprocess
import sys

import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def data_dir(exp, tag):
    return ensure_dir(os.path.join(PROJECT_ROOT, "data", exp, tag))


def results_dir(exp, tag):
    return ensure_dir(os.path.join(PROJECT_ROOT, "results", exp, tag))


def _git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=PROJECT_ROOT,
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def write_run_info(folder, args=None, extra=None):
    """写 run_info.json: 参数、代码版本、软件版本、机器、时间。每次运行都必须上传它。"""
    import scipy
    info = dict(time=datetime.datetime.now().isoformat(timespec="seconds"),
                command=" ".join(sys.argv), git_commit=_git_commit(),
                python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                platform=platform.platform(), machine=platform.node(),
                args=vars(args) if args is not None and hasattr(args, "__dict__") else args)
    if extra:
        info.update(extra)
    with open(os.path.join(folder, "run_info.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False, default=str)
    return info


def update_run_info(folder, **kw):
    path = os.path.join(folder, "run_info.json")
    info = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            info = json.load(f)
    info.update(kw)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False, default=str)


def save_solution(path, model, sol, extra=None, single_precision=False):
    """
    保存一个定态解 (npz)。single_precision=True 用 float32 存场 (用于上传的代表性解, 体积减半);
    需要再次作为牛顿初值或做传播时请用双精度版本 (data/ 目录)。
    """
    dt = np.float32 if single_precision else np.float64
    g = model.grid
    q = model.to_grid(np.asarray(sol["w"], complex))
    n_grid = model.index_grid(sol["n_t"] if sol.get("n_t") is not None
                              else model.thermal_index(np.abs(sol["w"]) ** 2))
    meta = dict(sample=model.sample.params, sample_name=model.sample.name, b=sol["b"],
                U=sol["power"], converged=sol.get("converged", True))
    if extra:
        meta.update(extra)
    np.savez_compressed(path, x=g.x, y=g.y, h=g.h, opt_mask=model.opt_mask,
                        th_type=model.sample.th_type, w_re=q.real.astype(dt), w_im=q.imag.astype(dt),
                        n=n_grid.astype(dt), b=sol["b"], U=sol["power"],
                        meta=json.dumps(meta, default=str))


def load_solution(path):
    d = np.load(path, allow_pickle=False)
    out = {k: d[k] for k in d.files}
    out["q"] = out["w_re"].astype(float) + 1j * out["w_im"].astype(float)
    out["meta"] = json.loads(str(out["meta"]))
    out["b"] = float(out["b"])
    out["U"] = float(out["U"])
    return out


def solution_on_model(saved, model):
    """
    把保存的解放到 model 的光学点上。若网格相同直接取值, 否则做三次插值 (网格加密/样品略变时用,
    插值后必须再跑一次牛顿迭代)。
    """
    from .diagnostics import interp
    from .geometry import Grid
    same = (len(saved["x"]) == len(model.grid.x) and len(saved["y"]) == len(model.grid.y)
            and np.allclose(saved["x"], model.grid.x) and np.allclose(saved["y"], model.grid.y))
    if same:
        return saved["q"][model.opt_mask]
    g = Grid.__new__(Grid)
    g.x, g.y, g.h = saved["x"], saved["y"], float(saved["h"])
    return interp(saved["q"], g, model.xo, model.yo)


def write_csv(path, rows):
    """rows: dict 列表 (键取所有行的并集)。"""
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow(r)


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))
