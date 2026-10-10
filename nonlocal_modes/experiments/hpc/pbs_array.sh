#!/bin/bash
# ============================================================================
# pbs_array.sh —— PBS Pro / OpenPBS 作业数组版本 (集群用 qsub 而不是 sbatch 时使用)
#
# 在 nonlocal_modes/experiments 目录下, 按阶段依次提交 (等上一阶段全部结束再提交下一阶段):
#     qsub -J 1-63 -v LIST=hpc/jobs/stage1.txt hpc/pbs_array.sh
#     qsub -J 1-48 -v LIST=hpc/jobs/stage2.txt hpc/pbs_array.sh
#     qsub -J 1-4  -v LIST=hpc/jobs/stage3.txt hpc/pbs_array.sh
# (任务数 = 清单文件的行数: wc -l hpc/jobs/stage1.txt)
# Torque (老式 PBS): 把 "-J" 换成 "-t", 并把下面的 PBS_ARRAY_INDEX 换成 PBS_ARRAYID。
# ============================================================================
#PBS -N nlmodes
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=48:00:00
#PBS -j oe

cd "$PBS_O_WORKDIR"
mkdir -p hpc/logs
# ---- 加载 Python 环境 (按实际情况修改) ----
# module load anaconda3
# source activate nl
IDX=${PBS_ARRAY_INDEX:?not running inside a PBS job array}
CMD=$(sed -n "${IDX}p" "$LIST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
STAGE=$(basename "$LIST" .txt)
START=$(date +%s)
eval "$CMD" > "hpc/logs/pbs_${STAGE}_${IDX}.out" 2>&1
RC=$?
END=$(date +%s)
echo -e "$IDX\t$RC\t$(( (END - START) / 60 ))min\t$CMD" >> "hpc/logs/status_${STAGE}.tsv"
exit $RC
