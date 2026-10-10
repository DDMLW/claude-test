#!/bin/bash
# ============================================================================
# submit_slurm.sh —— 一次性提交全部三个阶段 (阶段之间自动设置依赖)
#
# 在 nonlocal_modes/experiments 目录下:
#     python hpc/make_jobs.py            # 先生成任务清单 (默认 hq 精度)
#     bash hpc/submit_slurm.sh
# 可选环境变量:
#     SBATCH_EXTRA="-p cpu -A myaccount"   额外的 sbatch 参数 (分区、账户等)
#     MAXPAR=40                            同时运行的最大任务数 (默认不限制)
#     SBATCH=yhbatch                       天河等使用 yhbatch 的系统
# 例:  SBATCH_EXTRA="-p normal" MAXPAR=50 bash hpc/submit_slurm.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=${SBATCH:-sbatch}
EXTRA=${SBATCH_EXTRA:-}
THROTTLE=""
if [ -n "${MAXPAR:-}" ]; then THROTTLE="%${MAXPAR}"; fi
mkdir -p hpc/logs

submit () {   # submit <stage-file> [dependency-job-id]
    local list=$1 dep=${2:-} n jid
    n=$(grep -c . "$list")
    if [ -n "$dep" ]; then
        jid=$($SB --parsable $EXTRA --dependency=afterany:"$dep" --array=1-"$n""$THROTTLE" hpc/slurm_array.sh "$list")
    else
        jid=$($SB --parsable $EXTRA --array=1-"$n""$THROTTLE" hpc/slurm_array.sh "$list")
    fi
    echo "${jid%%;*}"
}

for k in 1 2 3; do
    [ -f "hpc/jobs/stage$k.txt" ] || { echo "hpc/jobs/stage$k.txt missing: run python hpc/make_jobs.py first"; exit 1; }
done
J1=$(submit hpc/jobs/stage1.txt)
J2=$(submit hpc/jobs/stage2.txt "$J1")
J3=$(submit hpc/jobs/stage3.txt "$J2")
echo "submitted: stage1=$J1  stage2=$J2 (after stage1)  stage3=$J3 (after stage2)"
echo "$(date '+%F %T') stage1=$J1 stage2=$J2 stage3=$J3" >> hpc/logs/submitted.txt
echo "monitor:  squeue -u \$USER      progress:  python hpc/status.py"
