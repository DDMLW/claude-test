"""
run_local.py —— 没有作业调度系统 (SLURM/PBS) 时, 在一台多核服务器上并行执行任务清单
====================================================================================
按依赖关系调度四个清单 (每个任务单核):
    stage1 与 indep 同时开始;  stage1 全部结束后开始 stage2;  stage2 与 indep 都结束后运行 final。
每个任务的屏幕输出写入 hpc/logs/local_<清单>_<行号>.out, 结束状态追加到 hpc/logs/status_<清单>.tsv
(与 SLURM 版本格式相同, 可用 hpc/status.py 查看)。

用法 (在 nonlocal_modes/experiments 目录下):
  python hpc/make_jobs.py
  nohup python hpc/run_local.py -n 32 > hpc/logs/run_local.out 2>&1 &     # 后台运行, 断开 ssh 也不停
  python hpc/run_local.py -n 8 --list stage2 --only 3 7     # 只重跑某个清单的第 3、7 行 (status.py 会给出命令)
  第二轮 (make_jobs.py --phase 2 之后):  nohup python hpc/run_local.py -n 32 --phase 2 > hpc/logs/run_local_p2.out 2>&1 &
  (先运行 p2_main 的全部任务, 结束后运行 p2_final)
注意: -n 不要超过物理核数; 每个任务最多约 6 GB 内存 (多数 1-2 GB), n x 内存不能超过机器内存。
"""

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, wait

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
LISTS = ("stage1", "indep", "stage2", "final", "p2_main", "p2_final")


def run_one(name, idx, cmd):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", MPLBACKEND="Agg", PYTHONUNBUFFERED="1")
    log = os.path.join(HERE, "logs", f"local_{name}_{idx}.out")
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%F %T')}] {name} task {idx}: {cmd}\n")
        f.flush()
        rc = subprocess.call(cmd, shell=True, cwd=EXP, env=env, stdout=f, stderr=subprocess.STDOUT)
    minutes = (time.time() - t0) / 60
    with open(os.path.join(HERE, "logs", f"status_{name}.tsv"), "a", encoding="utf-8") as f:
        f.write(f"{idx}\t{rc}\t{minutes:.0f}min\t{cmd}\n")
    print(f"[{time.strftime('%F %T')}] {name} #{idx} {'OK  ' if rc == 0 else 'FAIL'} "
          f"({minutes:.0f} min) {cmd}", flush=True)
    return rc


def read_list(name):
    path = os.path.join(HERE, "jobs", f"{name}.txt")
    if not os.path.exists(path):
        sys.exit(f"{path} not found: run python hpc/make_jobs.py first")
    with open(path, encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def submit(ex, name, only=None):
    cmds = read_list(name)
    todo = [(i, c) for i, c in enumerate(cmds, start=1) if only is None or i in only]
    print(f"===== {name}: {len(todo)} jobs queued =====", flush=True)
    return [(name, i, ex.submit(run_one, name, i, c)) for i, c in todo]


def report(futs):
    bad = [(name, i) for name, i, f in futs if f.result() != 0]
    for name in dict.fromkeys(n for n, _, _ in futs):
        nb = [i for n, i in bad if n == name]
        tot = sum(1 for n, _, _ in futs if n == name)
        print(f"{name} finished: {tot - len(nb)} OK, {len(nb)} failed {nb if nb else ''}", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-n", type=int, default=max(1, (os.cpu_count() or 2) - 1), help="同时运行的任务数")
    p.add_argument("--list", choices=LISTS, default=None, help="只运行这一个清单 (配合 --only 重跑失败的行)")
    p.add_argument("--only", nargs="+", type=int, default=None, help="只运行这些行号 (从 1 开始)")
    p.add_argument("--phase", type=int, choices=[1, 2], default=1, help="2: 第二轮清单 p2_main -> p2_final")
    a = p.parse_args()
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    with ThreadPoolExecutor(max_workers=a.n) as ex:
        if a.list:
            futs = submit(ex, a.list, a.only)
            wait([f for _, _, f in futs])
            report(futs)
        elif a.phase == 2:
            fm = submit(ex, "p2_main")
            wait([f for _, _, f in fm])
            report(fm)
            ff = submit(ex, "p2_final")
            wait([f for _, _, f in ff])
            report(ff)
        else:
            f1 = submit(ex, "stage1")
            fi = submit(ex, "indep")
            wait([f for _, _, f in f1])
            report(f1)
            f2 = submit(ex, "stage2")
            wait([f for _, _, f in f2 + fi])
            report(f2 + fi)
            ff = submit(ex, "final")
            wait([f for _, _, f in ff])
            report(ff)
    print("\nall done. check: python hpc/status.py")


if __name__ == "__main__":
    main()
