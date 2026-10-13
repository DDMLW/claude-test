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
python e2a_atlas.py --plot_only             # 拆分运行全部结束后: 从 data/ 中所有解族重建完整的 atlas_table.csv 与图谱
python e2a_atlas.py                         # 或者一次全跑 (约 1-3 小时)
python e2b_scaling.py                       # 约 1-2 小时
python e2b_scaling.py --h 0.15 --tag h015   # 加密网格复核标度律
```
> 拆分运行时，每次运行结束会用当时已有的解族重建 `atlas_table.csv` 和 atlas 图；全部拆分任务结束后请再运行一次 `--plot_only`，得到完整图谱。

**高精度（论文数据）**：在上面任一正式命令后加 `--hq`（结果写入 `hq/`），或在超算上用 [`../hpc/`](../hpc/README.md) 一次性拆分提交全部实验（推荐）。
hq 档：多边形/圆盘 h=0.1、矩形 h=0.15，b 步长 0.1，e2b 用 15 个 ρ（5.5–16，从 ρ=10 向两侧延拓；ρ 更大时 h=0.1 的网格需 8 GB 以上内存，且大 ρ 端 r_s 很小、用谐波拟合即可）。超算工具把 e2a 按（样品×m）拆成 36 个任务（每个 0.3–2 h），
并在粗一档网格（`--h_scale 1.5`，tag `hq_coarse`）上重算 8 个关键解族，最后由 `compare_tags.py` 生成
`hq/compare_hq_vs_hq_coarse.csv`：**中心荷、卫星数等整数量必须完全一致；U、ring_radius 的相对差应 < 1e-2**
（细网格结果的误差约为差值 × 0.8）。不一致的条目一般是卫星环 < 2h 的“未分辨”情形，请在 Progress.md 中列出。
> E3 需要读取 E2a 的原始数据，请保持相同的 `--tag`（默认 main），且不要删除 `nonlocal_modes/data/`。

## 输出与预期结果
### e2a_atlas.py
* `atlas_amplitude.png`、`atlas_phase.png`：行 = 样品（rect, poly3, poly4, poly5, poly6, disk），列 = m（1–6），b≈1 处。
  相位图中红圈 = +1 奇点，青圈 = −1，中心的叉号与数字 = 圆周绕数给出的中心荷。
* `atlas_table.csv`：每个 (样品, m) 一行，理论 vs 数值：`central_pred` / `central_num`、`satellites_pred` / `satellite_charge_num`、`mu_weight`（应 ≈1）、
  `r_s_fit`、`W_profile`（0.5h…3h 的绕数）、`resolved`（卫星环 ≥ 2h，结构可信）、`match`。
* `<样品>_m<m>_summary.csv`：沿 b 的完整解族（每个 b 一行）；`<样品>_m<m>_U_vs_b.png`；`<样品>_m<m>_b{0.5,1,2,4}.png`。
* 预实验（h=0.4，ρ=8，m=1–4）：已分辨的 22 个解族中 **21 个与理论一致**，包括矩形 m=4（中心 0 + 4 卫星，自动从正方形延拓得到）、
  五边形 m=4（中心 −1 + 5 卫星）、六边形 m=4（中心 −2 + 6 卫星）；三角形 m=2（r_s≈0.19）与五边形 m=3（r_s≈0.06）卫星环小于网格，标为 unresolved
  （h=0.2 时三角形 m=2 可分辨）；唯一不一致的是圆盘 m=4（网格诱导分裂，假象，见 E0d）。
* hq 云端测试（三角形 m=2，h=0.1 与 0.15 各 69 个解）：U 两套网格相对差 2.2e-3、ring_radius 3.5e-3、r_s_fit 2%；
  嬗变核随功率缩小（r_s_fit：b=0.2 时 0.9，b=2 时 0.09，b=6.8 时 0.017），所以网格读出的中心荷在 b ≲ 2 为 −1、高 b 端为 +2（未分辨，不是理论失效）；
  `compare_hq_vs_hq_coarse.csv` 中拓扑量的不一致全部出现在未分辨行（属正常），已分辨行 0 个不一致。
* **判读规则**（重要）：
  1. `resolved=False`（卫星环 < 2h）时，网格检测到的卫星位置不可信，只看 `r_s_fit` 与 `W_profile`；
  2. n=5,6 在 ρ=8 时物理各向异性与网格误差同量级（E0d），这两行只作参考，结论以 E5a 为准；
  3. 解族在 b 的下端停止（`low_end` 写着 “continuation stopped”）可能是**存在阈值**（理论预言 n | 2m 时出现），请记录阈值 b 值。

### e2b_scaling.py
* `scaling.png`：r_s（谐波拟合）vs ρ 的 log-log 图，图例给出拟合斜率与理论斜率；
* `scaling_fits.csv`：`slope_fit` vs `slope_pred`。预期：(3,2) −3，(4,3) −2，(3,3) −1，(4,4) −1。
  预实验（h=0.4）：(3,2) 拟合 **−3.009**；(3,3) 拟合 **−1.003**；(4,3) 拟合 −1.83（大 r_s 端偏离弱各向异性区，正式运行请看大 ρ 端的局部斜率）；
* `scaling_n<n>_m<m>.csv`：每个 ρ 的全部诊断量。

### e2c_export_fields.py（第五轮，论文 Fig.2 场图）
E2a 把每个解族的全部解以双精度存在 `data/E2_transmutation_atlas/<tag>/<样品>_m<m>_b<b>.npz`（不上传）。论文作图脚本只读 `results/`，
所以 E2c 对每个 (样品, m) 取 b 最接近 `--b` 的解，裁剪到光束区域（半宽 = 1.8 倍亮环半径 + 1）、降为单精度，写入
`results/E2_transmutation_atlas/<tag>/fields/<样品>_m<m>_b<b>.npz` 与 `fields/index.csv`。不做任何求解，几秒钟完成（hq 的 36×3 个场约 10 MB）。
```bash
python E2_transmutation_atlas/e2c_export_fields.py --hq --b 0.5 1.0 2.0 --job e2c_fields    # 超算上 (data/ 在那里)
python ../analysis/a9_atlas_fields.py --b 1.0                                                # 任意电脑: Fig.2 场图
```
quick 检查（2026-10-13，云端）：24 个场 0.2 MB，a9 画出的中心荷与卫星数逐格符合选择定则。
**如果超算上的 `data/E2_transmutation_atlas/hq/` 已被清理**，E2c 会对每个 (样品, m) 打印 “no saved solutions”：
此时用第一轮清单中 stage1 的 E2a 任务重跑需要的解族（每个 ≤ 1.4 h，结果与第一轮相同），再运行 E2c。

## 上传清单
`nonlocal_modes/results/E2_transmutation_atlas/main/` 或 `hq/`（以及 `hq_coarse`、`h015` 等复核目录）。Progress.md 中记录：
atlas_table 中 match 的个数与不一致的条目；各解族的 b 下限（阈值）；e2b 四个拟合斜率。
第五轮：`hq/fields/`（E2c 导出的场与 `index.csv`）与 `hq/log_e2c_fields.txt`。
