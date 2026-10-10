"""
run_local.py —— 没有作业调度系统 (SLURM/PBS) 时, 在一台多核服务器上并行执行任务清单
====================================================================================
依次执行 stage1 -> stage2 -> stage3; 每个阶段内同时运行 n 个任务 (每个任务单核)。
每个任务的屏幕输出写入 hpc/logs/local_<stage>_<序号>.out, 结束状态追加到 hpc/logs/status_<stage>.tsv
(与 SLURM 版本格式相同, 可用 hpc/status.py 查看)。

用法 (在 nonlocal_modes/experiments 目录下):
  python hpc/make_jobs.py
  nohup python hpc/run_local.py -n 32 > hpc/logs/run_local.out 2>&1 &     # 后台运行, 断开 ssh 也不停
  python hpc/run_local.py -n 32 --stages 2 3        # 只跑第 2、3 阶段
  python hpc/run_local.py -n 8 --only 3 7           # 只重跑 stage 清单中的第 3、7 行 (配合 --stages 使用)
注意: -n 不要超过物理核数; 每个任务最多约 3 GB 内存, 所以 n x 3 GB 不能超过机器内存。
"""

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)


def run_one(stage, idx, cmd):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", MPLBACKEND="Agg", PYTHONUNBUFFERED="1")
    log = os.path.join(HERE, "logs", f"local_{stage}_{idx}.out")
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%F %T')}] {stage} task {idx}: {cmd}\n")
        f.flush()
        rc = subprocess.call(cmd, shell=True, cwd=EXP, env=env, stdout=f, stderr=subprocess.STDOUT)
    minutes = (time.time() - t0) / 60
    with open(os.path.join(HERE, "logs", f"status_{stage}.tsv"), "a", encoding="utf-8") as f:
        f.write(f"{idx}\t{rc}\t{minutes:.0f}min\t{cmd}\n")
    print(f"[{time.strftime('%F %T')}] {stage} #{idx} {'OK  ' if rc == 0 else 'FAIL'} "
          f"({minutes:.0f} min) {cmd}", flush=True)
    return rc


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-n", type=int, default=max(1, (os.cpu_count() or 2) - 1), help="同时运行的任务数")
    p.add_argument("--stages", nargs="+", type=int, default=[1, 2, 3])
    p.add_argument("--only", nargs="+", type=int, default=None, help="只运行这些行号 (从 1 开始)")
    a = p.parse_args()
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    for k in a.stages:
        path = os.path.join(HERE, "jobs", f"stage{k}.txt")
        if not os.path.exists(path):
            sys.exit(f"{path} not found: run python hpc/make_jobs.py first")
        with open(path, encoding="utf-8") as f:
            cmds = [ln.strip() for ln in f if ln.strip()]
        todo = [(i, c) for i, c in enumerate(cmds, start=1) if a.only is None or i in a.only]
        print(f"\n===== stage{k}: {len(todo)} jobs, {a.n} in parallel =====", flush=True)
        with ThreadPoolExecutor(max_workers=a.n) as ex:
            rcs = list(ex.map(lambda ic: run_one(f"stage{k}", ic[0], ic[1]), todo))
        bad = [i for (i, _), rc in zip(todo, rcs) if rc != 0]
        print(f"stage{k} finished: {len(todo) - len(bad)} OK, {len(bad)} failed {bad if bad else ''}", flush=True)
    print("\nall done. check: python hpc/status.py")


if __name__ == "__main__":
    main()
