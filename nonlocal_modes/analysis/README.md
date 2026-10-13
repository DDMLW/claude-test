# analysis —— 由已上传结果表生成论文用的分析与图

这些脚本**只读取** `nonlocal_modes/results/` 中已上传的 CSV（不需要 `data/` 中的原始解，不做任何数值求解），
在任何电脑上几秒钟跑完，输出写到 `results/analysis/<tag>/`。师弟不需要运行它们；结果上传后由分析方运行并把结论写进 Progress.md。

图的统一样式在 `style.py`：固定顺序的色盲安全配色（蓝、橙、青、黄、品红、绿、紫、红），每个样品的颜色与标记固定
（矩形、三角形、正方形、五边形、六边形、圆盘），颜色从不是唯一的区分线索（同时用标记形状与图例），图中文字为英文。

| 脚本 | 读取 | 输出 | 回答的问题 |
|---|---|---|---|
| `a0_symmetry_transfer.py`（Fig.1） | E1 `e1a_index_anisotropy.csv`、`e1b_all.csv` | `symmetry_transfer.png`、`symmetry_transfer_summary.csv` | 对称性传递定律：dn_n·πn(n+1)/U 对各 n、ρ 塌缩到 (r/ρ)^n；系数与光束形状无关、对 U 线性（偏差 ≤ 1.6e-5）；保角半径相同的多边形与圆盘 U(b) 重合（m=0 ≤ 3.7e-4） |
| `a1_scaling_laws.py` | E2b `scaling_n*_m*.csv`、E5a `heating_n*_m*.csv` | `scaling_summary.csv`、`scaling_local_slopes.csv/png` | r_s ∝ ε^s 的指数是否等于牛顿多边形的预言？（渐近律，只在 r_s/r_ring → 0 时成立，所以给出局部斜率与渐近窗口 r_s/r_ring ≤ 0.2 内的拟合） |
| `a2_stability_map.py` | E3 `<样品>_m<m>_spectrum.csv`（及 `hq_coarse` 对照） | `stability_windows.csv`、`stability_map.png` | 各样品、各拓扑荷的失稳区间；低功率窗口的上端 b_c（Hamilton-Hopf 外推 g² ∝ b_c − b）；m ≥ 3 在高功率端的增长率 |
| `a3_e4_growth.py` | E4 `*_irreps.csv`、E2a 解族、E3 谱 | `e4_growth_vs_e3.csv/png` | E4 中 μ 权重漂移是不是“被种子放大的线性不稳定性”：非 μ 权重的指数增长率与同一定态（b(U)）的 E3 增长率比较 |
| `a4_window_mechanism.py`（Fig.4） | E3b 圆盘通道谱（`hq`、`hq_hib`）、各样品 m=2 谱、四档网格的三角形谱、E5c `hq_smallT` | `window_mechanism.png`、`hh_edge_check.csv` | 失稳窗口 = 圆盘 Bogoliubov 通道的 Krein 碰撞：(a) 低 b 段 L=−2（K+）与 L=+2/+4/+6（K−）的交叉 + 正方形/六边形失稳模频率；(b) 三角形 κ=0 对 L=−3/L=0 的连续追踪、失稳模频率 = 两支平均、Hamilton-Hopf 预言端点；(c) 各样品 m=2 增长率（矩形按保角半径换算到 ρ=8）；(d) 小 Θ 加热 γ/Θ 塌缩 |
| `a5_grid_convergence.py`（SI） | E3b 四档网格 `e3b_summary.csv` + a4 的 `hh_edge_check.csv` | `triangle_convergence.png/csv` | 三角形 m=2 窗口上端 b_c 与固定 b 处增长率随 h² 的外推；Hamilton-Hopf 端点 Δω_圆盘 = 2g(h) 与实测端点的对照；推出的耦合 g vs h² |
| `a6_dynamic_window.py`（Fig.5） | E4b 时间序列与汇总表（`hq`、`hq_h005`） | `dynamic_window.png`、`observable_window.csv` | (a) LG 发射后嬗变核 r_s(ξ)；(b) 非 μ 权重：稳定区代数增长（网格强迫）、定态发射与 1% 噪声不增长、失稳区 exp(2γξ)；(c) 可观测窗口 ξ_t50 → ξ_b |
| `a8_heated_window.py`（SI） | E5c `hq_hib`/`hq` 的 (3:2) 谱、半径 10 圆盘通道谱（`E3_stability/hq_disk10`）、`hq_hib_h0075/e5c_upper_edges.csv` | `heated_window.png/csv`、`heated_window_grid.csv` | C_3 加热 m=2 高 b 窗口 = 圆盘 κ=0 对 L=−3（K+）/L=0（K−）的 Hamilton-Hopf 碰撞：(a) 失稳模频率落在两支之间（≈ 平均 + O(Θ²)）；(b) c = sqrt(γ² + δ²/4)/Θ 对各 Θ 塌缩；(c) 预言的上端 δ(b_up) = 2cΘ（连续曲线）对 E5c 实测，加 h=0.075 的上端与 h² 外推 |
| `a7_twin_vortex.py`（Fig.6） | E5a `heating_n2_m2.csv`、E5c (2:2)、E5d `e5d_n2_m2_U15_summary.csv` 与 `fields/`、E5b `hq_omega/control_runs.csv` 与 `control_n2_*_tracks.csv`、E4b heat2 | `twin_vortex.png`、`twin_vortex_SI.png`、`twin_vortex_fields.png`、`twin_vortex_summary.csv` | 四极加热双涡旋：(a) r_s ∝ Θ^{1/2}；(b) (Θ, b) 稳定性 + U=15 工作点 + E5d 的 Θ_c；(c) 破坏 C_2 的平移模 ω_odd² 对 Θ 线性 → 叉式分岔 Θ_c=1.176；(d) 跟随相图（Ω, Θ）：E5b 传播的锁定 / 失锁点 + 一阶 / 二阶共振线 + Ω = −ω_odd(Θ)（破坏 C_2 的平移模在旋转系中零频），有 E5e 结果时叠加旋转坐标系同步定态的稳定段（同相支蓝、折返另一侧的支绿，上移 0.035）与失稳段（细线）——区间由 `_branch.csv` 按增长率阈值 `G_TH` = 2e-3（e 折长度 500）重建（`e5e_runs`），比汇总表的严格 1e-5 少算 6 个孤立的弱振荡泡；(e) r2_ratio 对零参数理论 ω_rot/(ω_rot − 2Ω)，叠加稳定同步定态的 r2（点线）；(f) 关断后自由共转 = ω_{−2}/2。`--e5e_tag` 可指向其他 tag 的 E5e 结果。SI：LG 发射中的形成、突然启动的滞后、间距；场图：各 Θ 的相位/振幅（`--fields_only` 只画场图） |
| `a9_atlas_fields.py`（Fig.2） | E2c 导出的 `E2_transmutation_atlas/<tag>/fields/*.npz`、E2a `atlas_table.csv` | `atlas_fields.png` | 嬗变图谱的场图：样品 × 拓扑荷，色相 = 相位、亮度 = sqrt 振幅，±1 奇点、中心荷、未分辨标记；圆盘对照 m ≥ 4 的中心分裂是方形网格自身的 C_4（标 “square-grid C_4”，Progress.md C3） |
| `fieldplot.py`（公用） | — | — | 场图的读取、HSV 相位-振幅着色（不用本身明暗变化的循环色图）、奇点标记、相位色条 |

运行（在 `nonlocal_modes/analysis/` 下）：
```bash
python a0_symmetry_transfer.py
python a1_scaling_laws.py --tag hq
python a2_stability_map.py --tag hq --coarse hq_coarse
python a3_e4_growth.py --tag hq
python a4_window_mechanism.py && python a5_grid_convergence.py     # a5 读取 a4 写出的 hh_edge_check.csv
python a6_dynamic_window.py
python a7_twin_vortex.py          # 有 E5d 场文件时同时画 twin_vortex_fields.png; --fields_only 只画场图
python a8_heated_window.py        # 需要第四轮的 hq_disk10 通道谱 (或 --channels 指定同格式文件)
python a9_atlas_fields.py --b 1.0 # 需要第五轮 E2c 导出的场 (results/E2_transmutation_atlas/hq/fields/)
```
第二轮实验（E3b、E4b、E5c）各自的 `--table_only` 已经生成汇总表与总图（`e3b_summary.csv`、`e4b_summary.csv`/`e4b_seed_scan.png`、
`e5c_summary.csv`/`e5c_windows.csv`/`e5c_stability_map.png`），结果上传后再决定是否需要新的分析脚本（例如把 E3b 的阈值与 E5c 的 Θ → 0 外推放在同一张图上）。

**2026-10 第一轮 hq 结果的分析结论**（详见 Progress.md C8、C9 与 §6.1）：
* a1：10 个指数（E2b 4 个、E5a 6 个）在渐近区全部与理论一致：最靠近渐近区的局部斜率与理论之差 ≤ 0.6%
  （例如 E2b (3,2) −3.0007 / 理论 −3，E2b (4,3) −1.988 / −2，E5a (3,2) 1.000002 / 1，E5a (4,3) 0.499 / 0.5）；全局拟合则因远离渐近区偏差可达 25%（E2b (4,4) −0.74）；
* a2：m=2 在所有样品中都有低功率失稳窗口（圆盘 b_c≈0.53，三角形 b_c≈2.7，五边形与圆盘相同），三角形 m=1 有窗口（b_c≈0.79），
  m ≥ 3 在所有样品、所有 b 都不稳定；
* a3：正方形与圆盘 m=2 中 μ 严格守恒（非 μ 权重 ≤ 1e-19，舍入误差级）；m ≥ 3 时 E4 的对称性破缺增长率不超过 E3 的线性增长率
  （E3 的增长率包含保持 C_n 对称的模式，是上界；例如圆盘 m=3 U=20：0.103 / 0.105）。但 m=2 的 LG 发射在定态稳定的 b 处也有慢增长
  （五边形/六边形 U=10：0.008/0.019；三角形 U=20：0.013，b(U)=2.80 恰在窗口边缘）——LG 光束不是定态，大幅呼吸时会扫过窗口内的 b；
  E4b 用宽度匹配的 LG、并把 U 取在远离窗口处来检验这一点。

**第二轮（E3b/E4b/E5c）的分析**（2026-10-12，结论见 Progress.md C8–C13）暂时直接由各实验的 `--table_only` 汇总表完成；
论文图阶段新增的 `a4`–`a7` 已在第三轮后实现（见上表与下文）。

**第三轮（2026-10-13）的分析结论**（详见 Progress.md C8–C13 与 §6.1“Claude 分析（第三轮）”）：
* a4/a5：三角形 m=2 宽窗口的机理定量成立——失稳模频率 = 圆盘 L=−3（K+）/L=0（K−）两支的平均（偏差 ≤ 0.14%），
  推出的耦合 g 在窗口内近似常数（h=0.05：0.017），边界判据 Δω_圆盘 = 2g 预言 b_c = 2.93，四档网格外推的实测上端 2.95–2.99；
  阈值随网格移动完全来自 g(h)（∝ h²）。小 Θ 加热（C_3，m=1）：g = 0.150 Θ（四个 Θ 偏差 ≤ 1.5%）。
* a6：稳定区的 LG 发射在 ξ≈2 内完成嬗变并保持到 600；失稳区的破缺距离由 exp(2γξ) 与种子决定（ξ_b≈70–280）。
* a7：双涡旋 r_s 渐近斜率 0.500；U=15 处 Θ≤1 全部稳定，低 b 实失稳窗口的上端随 Θ 线性上移（斜率 2.45），外推到 b(U=15) 给出 Θ_c≈1.2；
  关断后自由共转 0.1278 rad/ξ 对 L=−2 模预言 0.1262；慢转滞后 ≤ 0.05 rad；快转 Θ=1 在 ξ≈164 失锁。
* a8（云端圆盘通道谱，正式数据待第四轮 `hq_disk10`）：C_3 加热 m=2 高 b 窗口的耦合 c = sqrt(γ² + δ²/4)/Θ 对 Θ=1–3 塌缩到同一条 c(b)
  （b=2.8→5.8：0.019→0.009，∝ 亮环半径³ 精确到 1%，g ≈ 2.04 Θ (r_ring/R)³）；边界判据 δ_圆盘(b_up) = 2c(b_up)Θ 预言 3.358/4.207/4.920/5.545，实测 3.359/4.209/4.922/5.547。
* a0：传递定律塌缩到 (r/ρ)^n（n=3–6，ρ=6–18）；系数与光束形状、能流无关（偏差 ≤ 1.6e-5）；多边形与圆盘 U(b) 差 ≤ 3.7e-4。

**第四轮（2026-10-13）的分析结论**（详见 Progress.md C10、C13 与 §6.1“Claude 分析（第四轮）”）：
* a7（Fig.6 定稿）：E5d 正式数据 Θ_c = 1.176（ω_odd² = 0.0238 − 0.0202 Θ），E5c (2:2) Θ=1.2/1.3/1.5 的增长率落在同一直线上；
  Ω 扫描（52 个，用稳健的跟随指标重算）：锁定区间 Θ=0.25 −0.2…0.10、Θ=0.5 −0.2…0.06（0.03 失锁 = 二阶共振）、Θ=0.75 −0.1…0.01、Θ=1 −0.05…0.01；
  锁定时 r2_ratio 与 ω_rot/(ω_rot − 2Ω) 偏差 ≤ 4%（Θ ≤ 0.5、|Ω| ≤ 0.1、离开二阶共振）；h=0.075 关断共转 0.21904 对 h=0.1 的 0.2191。
* a8：正式 `hq_disk10` 通道谱与云端版本差 ≤ 4e-12；网格对照 Θ=1/2 的上端 3.359 → 3.431 → 外推 3.524、4.922 → 5.056 → 5.229（+5%/+6%）。
* a9 / a7 场图：读取第五轮导出的场（E2c、E5d `--fields_only`），quick 数据检查通过（选择定则的中心荷与卫星数逐格正确）。

**第五轮（2026-10-13）的分析结论**（详见 Progress.md C10 与 §6.1“Claude 分析（第五轮）”）：
* a7 (d)：E5e 9 个 Θ 的同步定态稳定区间（同相支含 Ω=0：Θ=0.125 −0.148…+0.029 → Θ=1.1 −0.043…+0.007）；反向边界与 −ω_odd 线在谱点间隔内重合
  （κ=1 平移模在旋转系中零频）；E5b 60 个点中 52 个与 E5e 一致，其余 8 个是 C_2 对称传播看不到的 κ=1、弱失稳或斜坡到不了的稳定态（第六轮加种子检验）。
* a9：Fig.2 用 108 个 hq 场重画（b≈1 的 36 格）；a7 场图用 6 个 hq 场重画。
