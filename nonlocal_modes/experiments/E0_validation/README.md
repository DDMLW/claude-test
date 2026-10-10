# E0 —— 代码验证（必须最先完成）

**目的**：证明代码可靠——与解析解一致、网格收敛、能复现导师已发表的两篇论文、并量化“网格假象”的底线。
审稿人一定会问这些问题，论文的 Methods / Supplement 直接用这里的结果。

## 运行顺序
```bash
cd nonlocal_modes
python -m pytest -q tests                  # 先跑单元测试, 必须全部通过 (约 5 s), 否则停下来报告
cd experiments/E0_validation
python e0a_convergence.py --quick          # 冒烟测试, 每个约 1-5 分钟
python e0b_twin_vortex.py --quick
python e0c_surface_solitons.py --quick
python e0d_grid_artifacts.py --quick
python e0a_convergence.py                  # 正式运行
python e0b_twin_vortex.py
python e0c_surface_solitons.py
python e0d_grid_artifacts.py
```
**高精度（论文数据）**：在上面任一正式命令后加 `--hq`（结果写入 `hq/`），或在超算上用 [`../hpc/`](../hpc/README.md) 一次性拆分提交全部实验（推荐）。
hq 档：e0a 网格加到 h=0.05（Poisson）/0.1（D_s、传递系数）、dz 扫到 0.00625；e0b/e0c 用 h=0.15、传播到 ξ=1000、
稳定性 shift 扫描 ω=0–2（9 个 shift）；e0d 网格加到 h=0.1。单核耗时：e0a 约 1 h（云端实测 40 分钟），e0b、e0c 各 4–10 h，e0d 约 1–2 h。
四个正式运行可以开四个终端同时跑。预计耗时（普通台式机）：e0a 15–30 min，e0b 1–2 h，e0c 1–2 h，e0d 0.5–1 h。

## 各脚本说明与预期结果

### e0a_convergence.py —— 收敛性
输出：`convergence_table.csv`、`convergence.png`
| 检验 | 预期 |
|---|---|
| [1] 圆盘高斯热源，与解析解比较 | 误差随 h 二阶下降（h 减半误差约 1/4）；h=0.1 时 ~1e-3 |
| [2] 对称性传递系数 数值/理论 | n=3,5,6：1.0000±1e-4；n=4 在 h≥0.3 时略偏（≤0.5%），h→0 趋于 1 |
| [3] 双涡旋 D_s(b=2)，20×40 矩形 | 随 h 收敛到约 **1.51**（h=0.4: 1.489, h=0.25: 1.508） |
| [4] 圆盘中 m=2 涡旋的 \|c₋₂/c₂\|（纯网格假象） | 九点格式比五点格式小 60–360 倍 |
| [5] 定态传播 ξ=20 的偏离度 | dz 减半，偏离度约降为 1/4；能流漂移 < 1e-12 |

### e0b_twin_vortex.py —— 复现 Ye et al., OL 35, 628 (2010)
输出：`e0b_fig1_m*_b*.png`（对应原文 Fig.1）、`e0b_fig1_table.csv`、`e0b_m2_20x40_summary.csv`、
`e0b_fig2a_fig3b.png`、`e0b_vs_deta.csv`、`e0b_fig2b_fig3a.png`、`e0b_prop_*`、`e0b_shape_perturbation_*`、`e0b_propagation.csv`

预期（**请导师对照原文图核对数值**）：
* m 个单位奇点沿**短轴 x（文中 η）**排成一行：m=2 中心无奇点 + 2 个卫星；m=3 中心 +1 + 2 个卫星；m=4 中心 0 + 4 个（m=4 若直接求解失败，脚本会自动从 40×40 正方形沿 d_η 延拓过来）；
* W_η ≥ W_ζ，d_η → 40（正方形）时二者相等；
* D_s 随 b 增大单调减小；随 d_η 减小而增大，在 20–40 之间变化最快，d_η < 15 后趋于饱和；正方形中 D_s → 0；
* 线性稳定性：m=1、2 稳定（max Re λ < 1e-6），m=3 不稳定（预实验增长率约 0.1）；
* 传播：m=1、2 偏离度保持在噪声水平（~0.1），m=3 在 ξ 几十到一两百内明显破坏；
* 20×40 的解放进 21×39 样品：两奇点连线角度周期振荡、不持续转动（`e0b_shape_perturbation_angle.png`）。

### e0c_surface_solitons.py —— 复现 Ye, Kartashov, Torner, PRA 2008
输出：`e0c_{dipole,vortex,bound_state}_b3.png`、`e0c_b3_states.csv`、`e0c_*_family_summary.csv`、`e0c_fig4ab.png`、`e0c_fig4c_phase.png`、`e0c_prop_*`、`e0c_propagation.csv`

预期：
* 偶极（节线垂直于界面）、表面涡旋（亮斑在 η 轴上、离界面远的更亮）、双涡旋束缚态（两个同号奇点沿 η 轴）都能在 b=3 求出；
  涡旋与束缚态的初值较难，脚本会自动用“自洽线性模式”初值搜索（会多花些时间，属正常）；
* U(b) 单调增，b→0 时 U→0；
* 表面涡旋的相位绕奇点非均匀增加（`e0c_fig4c_phase.png` 偏离对角线）；
* **注意（预实验发现，需要正式结果确认）**：表面偶极在 b=3 有一个很弱的振荡型线性不稳定（Re λ≈0.009，Im λ≈0.355，随网格加密略减小），
  但长距离传播（ξ=1200）中扰动非线性饱和，偶极保持结构、只有极点振幅的小幅振荡——与原文“稳定、两极振幅小幅振荡”的描述一致。
  请在 Progress.md 中如实记录正式运行的 max Re λ 和传播结果，供导师判断。

### e0d_grid_artifacts.py —— 网格假象的底线
输出：`e0d_grid_artifacts.csv`

* (a) 圆盘（物理上不应嬗变）：m=2、3 中心荷应保持 +2、+3；m=4、5 可能出现网格诱导的分裂，记录其随 h 的变化——这是“假象底线”。
  hq 云端测试：m=4 分裂半径 h=0.3/0.2/0.15/0.1 时为 0.45/0.35/0.23/0.28，**h=0.1 时并不更小**（m=4 是 r_s∝ε^{1/4} 的超敏感情形，
  圆形边界离散误差的 C_4 分量随 h 不单调），所以 m≥4 的任何结论都要与同 h 的圆盘对照比较；
* (b) 旋转检验：样品相对网格转 0°、7°、15°。**物理结构应随样品刚性转动**：中心荷相同，卫星半径相同，
  卫星相对样品的方位角（csv 中 `sat_angles_rel_mod`）相同。(3,2) 应通过；(5,3)、(6,4) 在 ρ=8 时预期**不**通过（物理各向异性太弱），这正是 n=5、6 改用 E5a 强加热来研究的原因。

## 上传清单
整个 `nonlocal_modes/results/E0_validation/main/`（或 `hq/`）目录，并在 Progress.md 运行记录中写明：
1. 单元测试是否全部通过；
2. e0a 表格中 [3] 的 D_s 收敛值；
3. e0b 中 m=1、2、3 的 max Re λ 与传播结论；d_η 扫描是否完整；
4. e0c 三个态是否都求出、表面偶极的 max Re λ 与传播结论；
5. e0d 旋转检验哪些通过、哪些不通过。
