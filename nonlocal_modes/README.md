# 非局域非线性材料中特殊模态的设计 —— 数值实验代码

**课题**：边界对称性主导的非局域涡旋孤子拓扑荷嬗变——从对称性传递定律到热边界可重构调控。
研究方案见 [`docs/research_plan.md`](docs/research_plan.md)，理论推导见 [`docs/theory.md`](docs/theory.md)，
数值方法与已知陷阱见 [`docs/numerics.md`](docs/numerics.md)，项目进度与运行记录见仓库根目录 [`Progress.md`](../Progress.md)。

所有实验都只需要电脑（Python），不需要线下实验。

---

## 1. 安装（一次即可）
需要 Python ≥ 3.9（推荐 Anaconda / Miniconda）。
```bash
git clone <仓库地址>
cd <仓库目录>
git checkout claude/nice-planck-5809ej        # 代码所在分支 (若导师已合并到 main, 用 main)
cd nonlocal_modes
pip install -r requirements.txt
python -m pytest -q tests                     # 24 个单元测试必须全部通过 (约 5 秒)
```
Windows 用户在 “Anaconda Prompt” 里执行同样的命令即可。

## 2. 跑实验（师弟的工作流程）
1. **先跑一遍所有 quick 模式**（粗网格冒烟测试，约 1 小时，确认在你的电脑上都能运行）：
   ```bash
   cd experiments
   python run_all.py --quick
   ```
   最后会打印每个脚本 OK / FAIL。quick 结果不用上传。
2. **按顺序做正式运行**：E0 → E1 → E2 → E3 → E4 → E5。每个实验目录里的 `README.md` 写明了
   运行命令、耗时、输出文件、**预期结果**和**上传清单**。E3 依赖 E2a 的数据（同一台电脑、同一 tag）。
   长时间任务可以开多个终端并行（脚本都支持 `--samples`、`--ms` 等参数拆分）。
3. **上传结果**：见 [`results/README.md`](results/README.md)。上传 `results/<实验>/<tag>/`，并在 `Progress.md` 的运行记录表里加一行。
4. 有任何报错或与预期不符：不要改代码，上传 `log.txt` 并在 Progress.md 里写清楚。

### 论文数据：hq 高精度档 + 超算批量运行
所有脚本都支持 `--hq`（h=0.1、矩形/界面样品 0.15，牛顿判据 1e-12，dz=0.01，b 步长减半，更长传播、更密稳定性扫描；结果在 `hq/` 下）。
计算量约为普通档的 10–30 倍，**请在超算上用 [`experiments/hpc/`](experiments/hpc/README.md) 的工具一次性拆分提交**（约 115 个单核任务，三个阶段，
SLURM/PBS/无调度服务器均可），其中自动包含粗一档网格的收敛对照（`hq_coarse`）。为什么不把收敛判据设成 1e-20：见 hpc/README.md 第 0 节与 docs/numerics.md 第 8 节。

| 实验 | 目录 | 内容 | 正式运行耗时(估, main 档) |
|---|---|---|---|
| E0 | `experiments/E0_validation` | 代码验证：解析解、收敛性、复现导师两篇论文、网格假象底线 | 3–5 h |
| E1 | `experiments/E1_symmetry_transfer` | 对称性传递定律、保角半径普适性 | 0.5–1 h |
| E2 | `experiments/E2_transmutation_atlas` | 拓扑荷嬗变“周期表”与标度律（核心） | 2–5 h |
| E3 | `experiments/E3_stability` | 线性稳定性谱 + 加噪声传播 | 4–10 h |
| E4 | `experiments/E4_dynamic_transmutation` | LG 光束动态嬗变（实验可验证的预言） | 1–3 h |
| E5 | `experiments/E5_thermal_reconfiguration` | 边界温度可重构调控 | 1–3 h |

## 3. 代码结构
```
nonlocal_modes/
├── nlmodes/                 核心数值包 (已通过单元测试)
│   ├── geometry.py          样品几何: 矩形/正多边形/圆盘/椭圆/绝热界面, 边界温度图样
│   ├── operators.py         各向同性九点差分 + Gibou 对称二阶边界的热传导算子
│   ├── model.py             非局域模型 (折射率 = 热传导方程的解)
│   ├── newton.py            定态牛顿迭代 (复/实, 定 b/定功率) + 参数延拓
│   ├── stability.py         线性稳定性 (稀疏 shift-invert Arnoldi; 稠密法)
│   ├── propagate.py         传播 (Crank-Nicolson + Strang 分步, 严格保能流)
│   ├── diagnostics.py       相位奇点、绕数剖面、角谐波、C_n 不可约表示投影
│   ├── guesses.py           初值 (LG、多极、自洽线性模式、理论指导的嬗变初值)
│   ├── theory.py            理论公式 (对称性传递系数、选择定则、牛顿多边形)
│   ├── io_utils.py          数据保存/读取、run_info.json
│   └── plotting.py          标准作图
├── experiments/             E0-E5 实验脚本, 每个目录有 README.md; run_all.py 一键运行
│   ├── common.py            共用工具 + 精度档位 (quick / main / hq)
│   ├── compare_tags.py      两个 tag 的结果表逐列比较 (网格收敛性检验)
│   └── hpc/                 超算批量运行: make_jobs.py, submit_slurm.sh, pbs_array.sh, run_local.py, status.py
├── tests/test_core.py       单元测试
├── results/                 上传的结果 (每个实验一个目录, 正式运行在 main/ 或 hq/ 下)
├── data/                    原始数据 (本地保存, 不上传)
└── docs/                    研究方案、理论、数值方法、参考文献
```

## 4. 模型与导师原始 MATLAB 代码的关系
导师提供的 `Relax2Dlocal` 是**局域**非线性 + 外势的牛顿迭代（雅可比只含对角非线性项）。本项目把同样的牛顿思路推广到
**非局域热非线性**：把折射率 n 也当作未知量，雅可比矩阵仍是稀疏的（3×3 分块），用稀疏 LU 求解；方程与归一化完全沿用
导师 PRA 2008 / OL 2010 两篇论文。为了研究微弱的 C_n 各向异性，又加入了各向同性九点格式、真实边界的二阶处理等改进（见 docs/numerics.md）。
用 Python 而不是 MATLAB 的原因：免费、可在任意电脑运行，且所有代码都已在云端逐一测试。
