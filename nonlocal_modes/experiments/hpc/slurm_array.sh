#!/bin/bash
# ============================================================================
# slurm_array.sh —— SLURM 作业数组: 第 i 个数组任务执行任务清单文件的第 i 行
#
# 一般不需要手动调用, 由 submit_slurm.sh 提交。手动提交 (在 nonlocal_modes/experiments 目录下):
#     sbatch --array=1-44 hpc/slurm_array.sh hpc/jobs/stage1.txt
# 只重跑失败的第 3、7 个任务:
#     sbatch --array=3,7 hpc/slurm_array.sh hpc/jobs/stage1.txt
#
# 资源: 每个任务单核 (SuperLU 稀疏分解是单线程的, 多给核没有用), 内存 8G 足够
# (实测 hq 档最大的任务约 1.2-3 GB, 见 hpc/README.md)。时间上限按你们集群分区的最大值调整。
# ============================================================================
#SBATCH --job-name=nlmodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=hpc/logs/slurm_%A_%a.out

# ---- 在这里加载你们集群上的 Python 环境 (按实际情况修改/取消注释) ----
# module load anaconda3
# source activate nl
# --------------------------------------------------------------------------

LIST="$1"
if [ -z "$LIST" ] || [ ! -f "$LIST" ]; then
    echo "usage: sbatch --array=1-N hpc/slurm_array.sh hpc/jobs/<stage1|stage2|indep|final>.txt"; exit 2
fi
IDX="${SLURM_ARRAY_TASK_ID:?not running inside a SLURM job array}"
CMD=$(sed -n "${IDX}p" "$LIST")
# 线性代数库只用单线程 (多个任务共享节点时互不抢核)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONUNBUFFERED=1

STAGE=$(basename "$LIST" .txt)
echo "[$(date '+%F %T')] $STAGE task $IDX on $(hostname): $CMD"
START=$(date +%s)
eval "$CMD"
RC=$?
END=$(date +%s)
echo "[$(date '+%F %T')] exit code $RC after $(( (END - START) / 60 )) min"
# 每个任务结束时追加一行状态, hpc/status.py 据此列出失败/未完成的任务
echo -e "$IDX\t$RC\t$(( (END - START) / 60 ))min\t$CMD" >> "hpc/logs/status_${STAGE}.tsv"
exit $RC
