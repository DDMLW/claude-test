#!/bin/bash
# ============================================================================
# submit_slurm.sh —— 一次性提交全部四个任务清单 (自动设置依赖关系)
#
#     stage1 (E2a) ──> stage2 (E3) ──┐
#     indep  (其余独立实验) ──────────┴──> final (汇总表、网格收敛对比)
#
# 在 nonlocal_modes/experiments 目录下:
#     python hpc/make_jobs.py            # 先生成任务清单 (默认 hq 精度)
#     bash hpc/submit_slurm.sh
# 可选环境变量:
#     SBATCH_EXTRA="-p cpu -A myaccount"   额外的 sbatch 参数 (分区、账户等)
#     MAXPAR=40                            每个作业数组同时运行的最大任务数 (默认不限制)
#     SBATCH=yhbatch                       天河等使用 yhbatch 的系统
# 例:  SBATCH_EXTRA="-p normal" MAXPAR=50 bash hpc/submit_slurm.sh
#
# 第二轮 (python hpc/make_jobs.py --phase 2 之后):
#     bash hpc/submit_slurm.sh p2          p2_main (全部独立) ──> p2_final (汇总表)
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SB=${SBATCH:-sbatch}
EXTRA=${SBATCH_EXTRA:-}
THROTTLE=""
if [ -n "${MAXPAR:-}" ]; then THROTTLE="%${MAXPAR}"; fi
mkdir -p hpc/logs

submit () {   # submit <list-file> [dependency-spec, 例如 afterany:123:456]
    local list=$1 dep=${2:-} n jid
    n=$(grep -c . "$list")
    if [ -n "$dep" ]; then
        jid=$($SB --parsable $EXTRA --dependency="$dep" --array=1-"$n""$THROTTLE" hpc/slurm_array.sh "$list")
    else
        jid=$($SB --parsable $EXTRA --array=1-"$n""$THROTTLE" hpc/slurm_array.sh "$list")
    fi
    echo "${jid%%;*}"
}

if [ "${1:-}" = "p2" ]; then
    for k in p2_main p2_final; do
        [ -f "hpc/jobs/$k.txt" ] || { echo "hpc/jobs/$k.txt missing: run python hpc/make_jobs.py --phase 2 first"; exit 1; }
    done
    JM=$(submit hpc/jobs/p2_main.txt)
    JF=$(submit hpc/jobs/p2_final.txt "afterany:$JM")
    echo "submitted: p2_main=$JM  p2_final=$JF (after p2_main)"
    echo "$(date '+%F %T') p2_main=$JM p2_final=$JF" >> hpc/logs/submitted.txt
    echo "monitor:  squeue -u \$USER      progress:  python hpc/status.py"
    exit 0
fi

for k in stage1 stage2 indep final; do
    [ -f "hpc/jobs/$k.txt" ] || { echo "hpc/jobs/$k.txt missing: run python hpc/make_jobs.py first"; exit 1; }
done
J1=$(submit hpc/jobs/stage1.txt)
JI=$(submit hpc/jobs/indep.txt)
J2=$(submit hpc/jobs/stage2.txt "afterany:$J1")
JF=$(submit hpc/jobs/final.txt "afterany:$J2:$JI")
echo "submitted: stage1=$J1  indep=$JI  stage2=$J2 (after stage1)  final=$JF (after stage2 and indep)"
echo "$(date '+%F %T') stage1=$J1 indep=$JI stage2=$J2 final=$JF" >> hpc/logs/submitted.txt
echo "monitor:  squeue -u \$USER      progress:  python hpc/status.py"
