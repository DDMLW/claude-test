"""
status.py —— 汇总超算任务的完成情况: 哪些成功、哪些失败、哪些还没有结束
=======================================================================
读取 hpc/logs/status_stage*.tsv (每个任务结束时由 slurm_array.sh / pbs_array.sh / run_local.py 追加一行),
与 hpc/jobs/stage*.txt 对照。同一任务重跑多次时以最后一次为准。

用法 (在 nonlocal_modes/experiments 目录下):  python hpc/status.py
失败的任务: 先看日志 (hpc/logs/ 下对应的 .out 文件, 以及 results/<实验>/<tag>/log_<job>.txt),
排除原因后只重跑失败的行, 例如 SLURM:  sbatch --array=3,7 hpc/slurm_array.sh hpc/jobs/stage1.txt
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    for k in (1, 2, 3):
        jobs = os.path.join(HERE, "jobs", f"stage{k}.txt")
        if not os.path.exists(jobs):
            continue
        with open(jobs, encoding="utf-8") as f:
            cmds = [ln.strip() for ln in f if ln.strip()]
        last = {}
        st = os.path.join(HERE, "logs", f"status_stage{k}.tsv")
        if os.path.exists(st):
            with open(st, encoding="utf-8") as f:
                for ln in f:
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) >= 3 and parts[0].isdigit():
                        last[int(parts[0])] = (int(parts[1]) if parts[1].lstrip("-").isdigit() else -1, parts[2])
        ok = [i for i in range(1, len(cmds) + 1) if last.get(i, (None,))[0] == 0]
        failed = [i for i in range(1, len(cmds) + 1) if i in last and last[i][0] != 0]
        pending = [i for i in range(1, len(cmds) + 1) if i not in last]
        print(f"stage{k}: {len(cmds)} jobs | OK {len(ok)} | FAILED {len(failed)} | not finished {len(pending)}")
        for i in failed:
            print(f"   FAILED #{i} (exit {last[i][0]}, {last[i][1]}): {cmds[i - 1]}")
        if failed:
            print(f"   resubmit (SLURM): sbatch --array={','.join(map(str, failed))} "
                  f"hpc/slurm_array.sh hpc/jobs/stage{k}.txt")
            print(f"   or locally:       python hpc/run_local.py --stages {k} --only {' '.join(map(str, failed))}")
        if pending and len(pending) <= 12:
            print(f"   not finished: {pending}")


if __name__ == "__main__":
    main()
