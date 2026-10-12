# analysis —— 由已上传结果表生成论文用的分析与图

这些脚本**只读取** `nonlocal_modes/results/` 中已上传的 CSV（不需要 `data/` 中的原始解，不做任何数值求解），
在任何电脑上几秒钟跑完，输出写到 `results/analysis/<tag>/`。师弟不需要运行它们；结果上传后由分析方运行并把结论写进 Progress.md。

图的统一样式在 `style.py`：固定顺序的色盲安全配色（蓝、橙、青、黄、品红、绿、紫、红），每个样品的颜色与标记固定
（矩形、三角形、正方形、五边形、六边形、圆盘），颜色从不是唯一的区分线索（同时用标记形状与图例），图中文字为英文。

| 脚本 | 读取 | 输出 | 回答的问题 |
|---|---|---|---|
| `a1_scaling_laws.py` | E2b `scaling_n*_m*.csv`、E5a `heating_n*_m*.csv` | `scaling_summary.csv`、`scaling_local_slopes.csv/png` | r_s ∝ ε^s 的指数是否等于牛顿多边形的预言？（渐近律，只在 r_s/r_ring → 0 时成立，所以给出局部斜率与渐近窗口 r_s/r_ring ≤ 0.2 内的拟合） |
| `a2_stability_map.py` | E3 `<样品>_m<m>_spectrum.csv`（及 `hq_coarse` 对照） | `stability_windows.csv`、`stability_map.png` | 各样品、各拓扑荷的失稳区间；低功率窗口的上端 b_c（Hamilton-Hopf 外推 g² ∝ b_c − b）；m ≥ 3 在高功率端的增长率 |
| `a3_e4_growth.py` | E4 `*_irreps.csv`、E2a 解族、E3 谱 | `e4_growth_vs_e3.csv/png` | E4 中 μ 权重漂移是不是“被种子放大的线性不稳定性”：非 μ 权重的指数增长率与同一定态（b(U)）的 E3 增长率比较 |

运行（在 `nonlocal_modes/analysis/` 下）：
```bash
python a1_scaling_laws.py --tag hq
python a2_stability_map.py --tag hq --coarse hq_coarse
python a3_e4_growth.py --tag hq
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
论文图阶段计划新增：`a4` 窗口机理图（圆盘通道谱 L、Krein 符号 + 各样品失稳窗口与失稳模频率，含 Hamilton-Hopf 边界判据表）、
`a5` 三角形 m=2 网格收敛图（阈值与增长率 vs h²）、`a6` E4b 可观测窗口图（ξ_t50 与 ξ_b）。
