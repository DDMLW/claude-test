# 超算批量运行（hq 高精度档，论文数据用）

这一目录的工具把全部实验拆成 115 个**单核任务**，写成四个任务清单，按依赖关系一次性提交到超算：

```
stage1 (E2a 解族, 44 个) ──> stage2 (E3 稳定性 + 图谱汇总, 46 个) ──┐
indep  (其余独立实验, 19 个, 提交后立即开始) ───────────────────────┴──> final (汇总表 + 网格收敛对比, 6 个)
```
做法是“多开任务并行”，而不是“一个任务用很多核”：稀疏 LU 分解（SuperLU）是单线程的，给一个任务多分核没有加速效果。

## 0. 精度档位说明（先读）
| 档位 | 命令行 | 网格 h（多边形/圆盘） | 网格 h（矩形 20×40、界面样品） | 牛顿收敛判据 | 传播步长 dz | b 步长 | 用途 |
|---|---|---|---|---|---|---|---|
| quick | `--quick` | 0.4 | 0.4 | 1e-9 | 0.02 | 粗 | 冒烟测试，不上传 |
| main | （默认） | 0.2 | 0.2 | 1e-9 | 0.02 | 0.2 | 普通正式运行 |
| **hq** | `--hq` | **0.1** | **0.15** | **1e-12** | **0.01** | **0.1** | **论文数据** |

hq 档还会延长传播距离（E0b/E0c/E3 到 ξ=1000，E4 到 400，E5b 到 600），加密稳定性 shift 扫描（E3：ω≤4、间隔 0.1、每个 shift 24 个本征值），
增加参数点（E0a 网格加到 h=0.05、E2b 15 个 ρ（5.5–16，从 ρ=10 向两侧延拓）、E5a 每种情形 19 个 Θ（公比 √2））。

**为什么收敛判据不设成 1e-20**：双精度浮点数的机器精度约 2.2e-16，差分算子把舍入误差放大 1/h² 倍，
残差的舍入噪声下限约 1e-13（实测 h=0.1 时 4e-13）。判据设得比这更小，所有求解都会被判为“不收敛”。
论文结果的误差主要来自网格（二阶格式，误差 ∝ h²，h=0.1 时约 1e-3～1e-4）和传播步长，而不是牛顿判据；
所以 hq 档把算力花在加密网格、减小步长、扩大参数扫描上，并用粗一档网格（tag `hq_coarse`）做收敛性对照，给出每个结果的误差估计。
程序在残差停在舍入下限（< 1e-10 且不再下降）时会判为收敛，并在日志中注明 `round-off floor reached`，这是正常的。

## 1. 准备（只做一次）
```bash
git clone https://github.com/DDMLW/claude-test.git
cd claude-test
git checkout claude/nice-planck-5809ej          # 若导师已合并到 main, 用 main
cd nonlocal_modes
# Python 环境 (按集群情况; 例如 module load anaconda3 后):
conda create -n nl python=3.11 numpy scipy matplotlib pytest -y
conda activate nl
python -m pytest -q tests                       # 29 passed (约 5 秒, 登录节点上跑即可)
```
SLURM 的 `sbatch` 默认继承提交时终端里的环境变量，所以**先 `conda activate nl` 再提交**一般就够了；
如果计算节点上找不到 python 或 numpy，再把 `experiments/hpc/slurm_array.sh`（或 `pbs_array.sh`）开头“加载 Python 环境”的两行
改成你们集群的实际写法（例如 `module load anaconda3` 和 `source activate nl`）。PBS 默认不继承环境，需要改这两行。

建议先提交一个 quick 冒烟测试，确认计算节点上的环境没问题（约 40 分钟）：
```bash
cd experiments
mkdir -p hpc/logs
sbatch -n 1 --mem=4G --time=03:00:00 -o hpc/logs/quick_%j.out --wrap "python run_all.py --quick"
```
输出最后应全部是 `OK`。

## 2. 生成任务清单并提交
```bash
cd nonlocal_modes/experiments
python hpc/make_jobs.py                          # 默认 hq; 生成 hpc/jobs/ 下 stage1 (44), stage2 (46), indep (19), final (6) 四个清单
# 若计算节点上 python 不是 conda 环境里的那个, 用 --python 指定绝对路径:
#   python hpc/make_jobs.py --python ~/miniconda3/envs/nl/bin/python
```
**SLURM**（`sbatch`/`squeue`）：
```bash
bash hpc/submit_slurm.sh                         # 四个清单一次提交, 自动设置依赖
# 需要指定分区/账户, 或限制同时运行的任务数时:
SBATCH_EXTRA="-p cpu -A 你的账户" MAXPAR=60 bash hpc/submit_slurm.sh
# 天河等使用 yhbatch 的系统:
SBATCH=yhbatch bash hpc/submit_slurm.sh
```
**PBS**（`qsub`/`qstat`）：没有自动依赖，按顺序提交（`-J` 后的数字是清单行数，用 `wc -l hpc/jobs/*.txt` 核对）：
```bash
qsub -J 1-44 -v LIST=hpc/jobs/stage1.txt hpc/pbs_array.sh      # 立即提交
qsub -J 1-19 -v LIST=hpc/jobs/indep.txt  hpc/pbs_array.sh      # 立即提交
qsub -J 1-46 -v LIST=hpc/jobs/stage2.txt hpc/pbs_array.sh      # stage1 全部结束后 (python hpc/status.py 查看)
qsub -J 1-6  -v LIST=hpc/jobs/final.txt  hpc/pbs_array.sh      # stage2 与 indep 都结束后
```
**没有作业调度系统的多核服务器**：
```bash
nohup python hpc/run_local.py -n 32 > hpc/logs/run_local.out 2>&1 &     # -n 不超过物理核数, 且 n x 3GB < 内存
```

## 3. 四个清单分别做什么
| 清单 | 任务 | 依赖 |
|---|---|---|
| stage1 | E2a 按（样品×m）拆成 36 个；网格收敛对照：8 个关键解族在粗一档网格上重算（tag `hq_coarse`） | 无 |
| stage2 | E2a 汇总图谱（`--plot_only`）；E3 按（样品×m）拆成 36 个（读取 stage1 的解族）；粗网格对照的 E2a 汇总与 8 个 E3 | stage1 全部结束 |
| indep | E0a–E0d、E1a、E1b、E2b、E5a、E5b 各 1 个；E4 按样品拆成 5 个（另 5 个加噪声，tag `hq_noisy`） | 无 |
| final | E3、E4 汇总表（`--table_only`）；`compare_tags.py` 比较 `hq` 与 `hq_coarse`（网格收敛性） | stage2 与 indep 全部结束 |

拆分出来的任务各自把日志写到 `results/<实验>/<tag>/log_<任务名>.txt`，互不覆盖。

## 4. 资源与耗时（单核；云端实测 + 估算）
云端实测（4 核服务器上 4 个任务同时运行，hq 网格）：
* E2a 一个解族（三角形 m=2，b=0.2–7.0 共 69 个解）：h=0.1 用时 11.5 分钟，粗网格 h=0.15 用时 3.4 分钟；残差全部 ≤ 4e-13；
* E3 一个 b 点的完整谱（41 个 shift × 24 个本征值）：h=0.1 约 8 分钟，h=0.15 约 3 分钟；
* 传播：每步 25–33 毫秒（ξ=1000 即 10 万步约 0.8–0.9 小时；E4 含奇点追踪约 30 毫秒/步）；
* 单个进程内存 1.2–1.3 GB（多边形 h=0.1、矩形 h=0.15、R=10 圆盘 h=0.1）。

| 任务 | 每个任务耗时（估） | 内存 |
|---|---|---|
| E2a（每个样品×m） | 0.2–1 h（矩形 m≥4 需从 40×40 正方形延拓，可能到 2 h） | ≤ 4 GB（正方形 7 万格点） |
| E3（每个样品×m） | 5–8 h（约 35 个 b 点 × 8 分钟 + 2 次 ξ=1000 传播） | ≤ 2 GB |
| E0a | 约 1 h（云端实测 40 分钟） | ≤ 6 GB（矩形 h=0.1 有 8 万格点） |
| E0b、E0c | 4–10 h（ξ=1000 传播 3–4 次；E0c 还有表面束缚态搜索） | ≤ 4 GB |
| E4（每个样品，6 个 run 到 ξ=400） | 约 2 h | ≤ 2 GB |
| E1b、E5a、E5b | 1–5 h | ≤ 2 GB |
| E2b | 3–8 h（ρ=16 时每次牛顿迭代约 30 s） | ≤ 6 GB（ρ=16 实测峰值 4.9 GB） |

`slurm_array.sh` 默认申请 `--mem=8G --time=48:00:00`（单个任务最长估计约 10 h，留有余量；若分区上限是 24 h 也够用）。总计约 300–400 核时；若能同时运行 60 个以上任务，
墙钟时间约 10–12 小时（最长的一条链：stage1 约 2 h + stage2 约 8–10 h；indep 中最长的 E0b/E0c 约 10 h 同时进行；final 几分钟）。实际耗时会写在 `hpc/logs/status_stage*.tsv` 中，请记录到 Progress.md。

## 5. 查看进度、处理失败
```bash
squeue -u $USER                  # SLURM 队列 (PBS: qstat -u $USER)
python hpc/status.py             # 每个清单: 成功 / 失败 / 未结束, 并打印失败任务的重提交命令
```
* 某任务失败：看 `hpc/logs/slurm_<作业号>_<序号>.out` 和 `results/<实验>/<tag>/log_<任务名>.txt` 的最后几十行；
  内存不足（`OOM`/`oom-kill`）就把 `slurm_array.sh` 里的 `--mem` 改成 16G；超时就把 `--time` 加大。然后只重提交失败的行，例如
  `sbatch --array=3,7 hpc/slurm_array.sh hpc/jobs/stage1.txt`。
* 某个 E2a 任务失败后重跑成功：还要重跑对应的 E3 任务（stage2 中同名的那一行），以及 stage2 第 1 行（图谱汇总）和 final 第 1 行（E3 汇总）。
* 不要改代码。解决不了的报错：上传日志，在 Progress.md 中写清楚命令和报错最后几行。

## 6. 跑完之后上传什么
```
nonlocal_modes/results/<每个实验>/hq/          全部上传 (论文数据)
nonlocal_modes/results/E2_transmutation_atlas/hq_coarse/   上传 (网格收敛对照)
nonlocal_modes/results/E3_stability/hq_coarse/             上传
nonlocal_modes/results/E4_dynamic_transmutation/hq_noisy/  上传
```
`data/` 目录不上传（很大），但请保留在超算上，之后可能还要读取。
重点结果都在：`results/E2_transmutation_atlas/hq/compare_hq_vs_hq_coarse.csv` 与 `results/E3_stability/hq/compare_hq_vs_hq_coarse.csv`
（网格收敛性：输出中 “resolved rows” 的拓扑量不一致数必须为 0 —— “unresolved rows” 是卫星环小于网格的行，不一致属正常；
连续量的相对差应 < 1e-2；E3 的 `unstable` 列不一致数应为 0），以及 `python hpc/status.py` 的输出 —— 请把它贴进 Progress.md 的运行记录。
上传方法见 `nonlocal_modes/results/README.md`。

## 7. 第二轮（phase 2，2026-10）：补充实验
第一轮（115 个任务）的分析结论与第二轮的理由见 `Progress.md` 第 6.1 节与 `docs/research_plan.md` 第 6.1 节。
第二轮共 **55 个计算任务 + 8 个汇总任务**，全部互相独立（不读取第一轮的 `data/`），只需两个清单：
```
p2_main (55 个, 提交后立即全部开始) ──> p2_final (8 个汇总表, 几分钟)
```
**不会**覆盖第一轮的清单和 `hpc/logs/status_stage1.tsv` 等状态文件；`python hpc/status.py` 会同时列出两轮。

```bash
cd nonlocal_modes/experiments
git pull                                         # 先同步代码 (第二轮新增 E3b / E4b / E5c)
python -m pytest -q ../tests                     # 应为 29 passed
python hpc/make_jobs.py --phase 2                # 生成 hpc/jobs/p2_main.txt (55) 与 p2_final.txt (8)
# SLURM:
bash hpc/submit_slurm.sh p2                      # 可加 SBATCH_EXTRA="-p ..." MAXPAR=60, 同第 2 节
# PBS:
qsub -J 1-55 -v LIST=hpc/jobs/p2_main.txt hpc/pbs_array.sh
qsub -J 1-8  -v LIST=hpc/jobs/p2_final.txt hpc/pbs_array.sh   # p2_main 全部结束后
# 无作业调度的多核服务器:
nohup python hpc/run_local.py -n 32 --phase 2 > hpc/logs/run_local_p2.out 2>&1 &
python hpc/status.py                             # 进度
```

| 实验 | 任务数 | 内容 | 单核耗时（估） |
|---|---|---|---|
| E3b | 15 | 8 个情形的失稳窗口（阈值、模式结构、非线性结局）；三角形 m=2 高 b 端在同一组 b 点上的网格收敛表 h=0.15/0.075/0.05（tag `hq_coarse`/`hq_h0075`/`hq_h005`，与 `hq` 一起按 h² 外推）与旋转 15° 检验；圆盘 m=2 阈值在 h=0.05 加密；圆盘 m=1、2 的通道分解谱 | 2–8 h；h=0.05 的两个约 7–11 h（**最先开始**） |
| E5c | 24 | 3 个情形（C_3 加热 m=2、C_4 加热 m=2、C_3 加热 m=1）× 8 个 Θ 的 (Θ, b) 稳定性 | 3–6 h |
| E4b | 16 | 14 个发射（稳定区、种子扫描、圆盘对照、加热样品）+ 2 个网格对照（`hq_coarse`、`hq_h005`） | 1–3 h；h=0.05 约 6 h |

总计约 250 核时；能同时运行 55 个任务时墙钟约 10–12 h（由 h=0.05 的 E3b 加密决定）。
内存：h=0.05 的任务峰值约 5 GB（云端实测 4.9 GB；`slurm_array.sh` 默认 8G 够用），其余 ≤ 2 GB。
（2026-10-12 修正了稳定性计算中相邻两个 shift 的 LU 同时占内存的问题，修正前 h=0.05 峰值 9.9 GB，会超过 8G。）

**跑完上传**（`git add nonlocal_modes/results Progress.md`，`data/` 与 `quick*` 自动忽略）：
```
results/E3_stability/hq/            新增 e3b_* 与 log_*         results/E3_stability/hq_coarse/, hq_h0075/, hq_h005/   新增 e3b_*
results/E5_thermal_reconfiguration/hq/   新增 e5c_*
results/E4_dynamic_transmutation/hq/     新增 e4b_*            results/E4_dynamic_transmutation/hq_coarse/, hq_h005/
```
并把 `python hpc/status.py` 的输出与总耗时贴进 Progress.md 第 6 节的运行记录。

## 8. 第三轮（phase 3，2026-10-12）：补漏 + 新实验
第二轮的分析与第三轮的理由见 `Progress.md` §6.1“Claude 分析（第二轮）”与 §8。共 **34 个计算任务 + 5 个汇总任务**，全部互相独立：
```bash
cd nonlocal_modes/experiments
git pull && python -m pytest -q ../tests        # 应为 29 passed
python hpc/make_jobs.py --phase 3                # 生成 hpc/jobs/p3_main.txt (34) 与 p3_final.txt (5)
bash hpc/submit_slurm.sh p3                      # PBS: qsub -J 1-34 -v LIST=hpc/jobs/p3_main.txt hpc/pbs_array.sh, 结束后 1-5 的 p3_final
python hpc/status.py
```
| 任务 | 数量 | 内容 | 单核耗时（估） |
|---|---|---|---|
| E5c C_3 加热高 b 端（`hq_hib`） | 5 | (3:2)，Θ=1–3，b=2.9–6.0 | 4–6 h |
| E5c 四极加热（`hq`） | 8 | (2:2)，Θ=0–1 | 3–7 h |
| E5c (3:1) 重跑（`hq`） | 8 | b 范围扩到 3.0（第二轮 Θ=3 无数据） | 3–7 h |
| E3b | 2 | rect:2 重跑（解族起步修正）；圆盘 m=2 通道谱 b=1.4–3.4（`hq_hib`） | 2–5 h |
| E5c 小 Θ 线性律（`hq_smallT`） | 5 | (3:1)，Θ=0–0.16，b=0.1–0.5、间隔 0.01 | 2–3 h |
| E5b 四极加热双涡旋 | 2 | Θ=0.5、1，关断/慢转/快转，ξ 到 600 | 1–2 h |
| E4b | 4 | 定态发射对照 2 个；四极加热中的 LG 发射 2 个 | 1–2 h |

**上传**：`results/E3_stability/{hq,hq_hib}`、`results/E5_thermal_reconfiguration/{hq,hq_hib,hq_smallT}`、`results/E4_dynamic_transmutation/hq`
中新增或更新的文件（(3:1) 的 `e5c_n3_m1_*` 会被重跑覆盖，这是预期的），并在 Progress.md §6 记录运行情况。

## 9. 第四轮（phase 4，2026-10-13）：双涡旋的旋转跟随（Fig.6）
第三轮的分析与第四轮的理由见 `Progress.md` §6.1“Claude 分析（第三轮）”与 §8、`docs/theory.md` §5.1。共 **61 个计算任务 + 4 个汇总任务**，全部互相独立：
```bash
cd nonlocal_modes/experiments
git pull && python -m pytest -q ../tests        # 应为 30 passed
python hpc/make_jobs.py --phase 4                # 生成 hpc/jobs/p4_main.txt (61) 与 p4_final.txt (4)
bash hpc/submit_slurm.sh p4                      # PBS: qsub -J 1-61 -v LIST=hpc/jobs/p4_main.txt hpc/pbs_array.sh, 结束后 1-4 的 p4_final
python hpc/status.py
```
| 任务 | 数量 | 内容 | 单核耗时（估） |
|---|---|---|---|
| E3b 圆盘通道谱（`hq_disk10`） | 1 | 半径 10 圆盘（与 E5 同一网格）m=2，b=2.8–6.2、Im 0.55–1.3：C_3 加热高 b 窗口的 Hamilton-Hopf 检验所需的圆盘频差 | 3–5 h |
| E5c 网格对照（`hq_hib_h0075`） | 2 | C_3 加热 m=2 高 b 窗口上端（Θ=1：b=3.0–4.2；Θ=2：b=4.4–5.6）在 h=0.075 上重算：上端随网格移动多少（云端预算 Θ=1：3.36 → ≈3.43） | 2–4 h |
| E5d 双涡旋的低频谱（`hq`） | 1 | U=15，Θ=0–1.5：转动/分裂模频率 ω_{−2}(Θ) → 预言 Ω_c = ω_{−2}/2 | 0.5–1 h |
| E5c 四极加热（`hq`） | 4 | (2:2)，Θ=1.1、1.2、1.3、1.5：U=15 处双涡旋稳定的上限（云端谱给出 Θ_c=1.18：预言 1.1 稳定，1.2/1.3/1.5 失稳） | 3–7 h |
| E5b 网格对照（`hq_h0075`） | 1 | Θ=1 的关断与快转（Ω=0.05）在 h=0.075 上重复：失锁距离（h=0.1 时 ξ≈164）与共转角速度是否不变 | 2–3 h |
| E5b Ω 扫描（`hq_omega`） | 52 | Θ=0.25/0.5/0.75/1 × Ω=+0.01…+0.06（间隔 0.01）、0.08、0.10、0.12、0.14 与 −0.05/−0.10/−0.20；角速度在 ξ≤100 内线性升起，传播到 ξ=600 | 每个约 0.5 h |

**上传**：`results/E5_thermal_reconfiguration/{hq,hq_omega,hq_h0075,hq_hib_h0075}`、`results/E3_stability/hq_disk10` 中新增或更新的文件（`hq` 的 e5c 汇总表会被 p4_final 重建，属正常），
并在 Progress.md §6 记录运行情况。自查：`hq_omega/control_runs.csv` 中每个 Θ 的正向 Ω 是否在某个值以上出现 `xi_lock_lost`
（失锁），负向是否全部不失锁；`r2_ratio` 是否随 Ω 单调（正向 > 1、反向 < 1）。
