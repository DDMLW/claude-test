# E2 —— 拓扑荷嬗变“周期表”与标度律（论文 Fig.2、Fig.3，核心结果）

**物理问题**：电荷 m 的涡旋孤子在 C_n 对称的均匀热样品中，相位奇点如何排布？

**理论预言**（docs/theory.md 第 3 节；`python -c "from nlmodes import theory; print(theory.selection_table())"` 可打印全表）：
* 角赝动量 μ = m mod n 守恒；中心荷 v = |l| 最小的允许谐波；其余 m − v 个单位荷以 C_n 对称“卫星环”分布；
* 卫星环半径 r_s ∝ ε^s，ε ∝ U/ρ^n，所以**固定 U 时 r_s ∝ ρ^(−n s)**；
* 例：矩形 m=2 → 中心 0 + 2 卫星（导师的双涡旋）；三角形 m=2 → 中心 **−1** + 3 个 +1；正方形 m=3 → 中心 **−1** + 4 个 +1；三角形 m=3 → 0 + 3；正方形 m=2 → 中心 +2（μ=n/2 简并情形）。

## 运行
```bash
cd nonlocal_modes/experiments/E2_transmutation_atlas
python e2a_atlas.py --quick                 # 约 15-30 分钟 (粗网格 h=0.4, m<=4)
python e2b_scaling.py --quick               # 约 5 分钟
# 正式运行 (可开多个终端按样品拆分, 结果写入同一目录, 最后 atlas 图需全部样品都算完后再跑一次只画图的命令, 见下):
python e2a_atlas.py --samples rect
python e2a_atlas.py --samples poly3 poly4
python e2a_atlas.py --samples poly5 poly6 disk
python e2a_atlas.py                         # 或者一次全跑 (约 1-3 小时)
python e2b_scaling.py                       # 约 1-2 小时
python e2b_scaling.py --h 0.15 --tag h015   # 加密网格复核标度律
```
> 拆分运行时，每次运行会覆盖 `atlas_table.csv` 和 atlas 图（只含本次的样品）；最终论文图请用一次全样品运行生成。
> E3 需要读取 E2a 的原始数据，请保持相同的 `--tag`（默认 main），且不要删除 `nonlocal_modes/data/`。

## 输出与预期结果
### e2a_atlas.py
* `atlas_amplitude.png`、`atlas_phase.png`：行 = 样品（rect, poly3, poly4, poly5, poly6, disk），列 = m（1–6），b≈1 处。
  相位图中红圈 = +1 奇点，青圈 = −1，中心的叉号与数字 = 圆周绕数给出的中心荷。
* `atlas_table.csv`：每个 (样品, m) 一行，理论 vs 数值：`central_pred` / `central_num`、`satellites_pred` / `satellite_charge_num`、`mu_weight`（应 ≈1）、
  `r_s_fit`、`W_profile`（0.5h…3h 的绕数）、`resolved`（卫星环 ≥ 2h，结构可信）、`match`。
* `<样品>_m<m>_summary.csv`：沿 b 的完整解族（每个 b 一行）；`<样品>_m<m>_U_vs_b.png`；`<样品>_m<m>_b{0.5,1,2,4}.png`。
* 预实验（h=0.4，ρ=8）：矩形 m=1–3、三角形 m=1–4、正方形 m=1–4、圆盘 m=1–3 **全部与理论一致**；
  六边形 m=4 出现预言的 6 个卫星环（r≈1.2），但中心 −2 在亚网格尺度被网格劈裂；圆盘 m=4 出现网格诱导分裂（假象，见 E0d）。
* **判读规则**（重要）：
  1. `resolved=False`（卫星环 < 2h）时，网格检测到的卫星位置不可信，只看 `r_s_fit` 与 `W_profile`；
  2. n=5,6 在 ρ=8 时物理各向异性与网格误差同量级（E0d），这两行只作参考，结论以 E5a 为准；
  3. 解族在 b 的下端停止（`low_end` 写着 “continuation stopped”）可能是**存在阈值**（理论预言 n | 2m 时出现），请记录阈值 b 值。

### e2b_scaling.py
* `scaling.png`：r_s（谐波拟合）vs ρ 的 log-log 图，图例给出拟合斜率与理论斜率；
* `scaling_fits.csv`：`slope_fit` vs `slope_pred`。预期：(3,2) −3，(4,3) −2，(3,3) −1，(4,4) −1。
  预实验（h=0.4）：(3,2) 拟合 **−3.009**；(4,3) 拟合 −1.83（大 r_s 端偏离弱各向异性区，正式运行请看大 ρ 端的局部斜率）；
* `scaling_n<n>_m<m>.csv`：每个 ρ 的全部诊断量。

## 上传清单
`nonlocal_modes/results/E2_transmutation_atlas/main/`（以及 `h015` 等复核目录）。Progress.md 中记录：
atlas_table 中 match 的个数与不一致的条目；各解族的 b 下限（阈值）；e2b 四个拟合斜率。
