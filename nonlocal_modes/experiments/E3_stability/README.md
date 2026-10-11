# E3 —— 嬗变态的稳定性（论文 Fig.4，最可能出现“意外发现”的部分）

**前提**：先完成 E2a，且使用相同的 `--tag`（默认 main）。本脚本读取 `nonlocal_modes/data/E2_transmutation_atlas/<tag>/` 中的解族，所以必须在跑 E2a 的同一台电脑上运行。

**物理问题**：边界对称性造成的嬗变结构（中心反号奇点 + 卫星环）会不会改变高阶涡旋的稳定性？
已知：圆柱样品 m≤2 稳定（Kartashov 2007）；矩形 m=2 稳定、m≥3 不稳定（Ye 2010）。**多边形中的结果未知**——无论稳定与否都是新结果，请如实记录。

## 运行
```bash
cd nonlocal_modes/experiments/E3_stability
python e3_stability.py --quick                         # 需先跑 E2a --quick; 约 10-20 分钟
# 正式运行很耗时 (每个解要扫 ~20 个 shift), 建议按样品拆开, 多开终端并行:
python e3_stability.py --samples poly3
python e3_stability.py --samples poly4
python e3_stability.py --samples rect
python e3_stability.py --samples disk
python e3_stability.py --samples poly5 poly6           # 参考用 (见 E2 README 判读规则)
# 只想先看线性谱、暂不做传播:
python e3_stability.py --samples poly3 --no_prop
# 拆分运行全部结束后重建总表:
python e3_stability.py --table_only
```
**高精度（论文数据）**：在上面任一正式命令后加 `--hq`（结果写入 `hq/`），或在超算上用 [`../hpc/`](../hpc/README.md) 一次性拆分提交全部实验（推荐）。
hq 档（需先跑 E2a `--hq`）：shift 扫描 ω=0–4、间隔 0.1（41 个 shift）、每个 shift 24 个本征值，b 间隔 0.2（E2a 步长 0.1、stride 2），
传播到 ξ=1000（dz=0.01）。单核每个（样品×m）解族约 5–12 h；超算工具会把它拆成 36 个任务，并对 8 个关键解族在粗网格（`hq_coarse`）上重算，
由 `compare_tags.py` 给出稳定性判定是否随网格收敛（`hq/compare_hq_vs_hq_coarse.csv` 中 `unstable` 列不一致的行数应为 0）。
可调参数：`--stride`（每隔几个 b 算一次谱，默认 2）、`--omega_max`（shift 扫描的最大虚部，默认 3）、`--zmax`（传播距离，默认 500）。
预计耗时：多边形每个解族约 10–30 分钟（含传播），矩形每个解族约 1–2 小时。

## 输出与预期结果
* `<样品>_m<m>_spectrum.csv`：每个 b 的 `max_growth`（= max Re λ）、最不稳定本征值的实部/虚部、是否不稳定（> 1e-4）；
* `<样品>_m<m>_stability.png`：左 U(b)，蓝点稳定、红叉不稳定；右 max Re λ vs b（对数坐标，虚线为判据 1e-4）；
* `<样品>_m<m>_prop_b*_snapshots.png / _timeseries.png`：加噪声（方差 0.01）传播的快照与偏离度；
* `<样品>_m<m>_propagation.csv`：`survived`（后半程最大偏离 < 0.5）；
* `stability_table.csv`：总表。

**预实验（h=0.4，粗网格，每族只算 4–7 个 b，仅供参考）**：
* 所有多边形（n=3–6）中 m=2 嬗变态在 b ≳ 2 线性稳定，但在低功率 b=0.5 不稳定（max Re λ 0.01–0.04，五边形仅 2e-4）；
  圆盘 m=2 在全部 b 稳定。→ 可能是“边界对称性诱导的低功率失稳”，**请重点在高分辨率下确认，并加密 b 网格找出失稳阈值**；
* 所有 m=3、4 态（多边形与圆盘）不稳定，增长率 0.11–0.19；
* 矩形 m=1、2 稳定，m=3、4 不稳定（与导师 2010 年结论一致）；圆盘 m=1、2 稳定，m=3、4 不稳定（与 Kartashov 2007 一致）；
* quick 模式只传播到 ξ=40，不足以让增长率 ~0.1 的不稳定完全发展，正式运行（ξ=500）才能看到破坏。

**hq 云端定点测试（三角形 m=2，`--bs 0.5 1 3`）**：低功率失稳在加密网格下**确认存在且已收敛**——
b=0.5：max Re λ = 0.04110（h=0.1）/ 0.04104（h=0.15），Im λ = 0.118；b=1.0：0.0415 / 0.0413（Im 0.162）；b=3.0：稳定（< 1e-6）。
hq 正式结果中三角形 m=2 应在 b ≈ 1–3 之间出现稳定性转变（b 间隔 0.2 可定出阈值），这一族请重点看。

**判读规则**：
1. 线性谱判不稳定但传播中 `survived=True` 且偏离度有界振荡：记为“弱不稳定、非线性饱和”（不是数值错误，见 Progress.md 中表面偶极的例子）；
2. 线性谱判稳定但传播破坏：极少见，可能是 shift 扫描漏掉了本征值，请用 `--omega_max 5` 重算该族并报告；
3. 若发现某个 m ≥ 3 的嬗变态在一段 b 内**稳定**，这是重要结果：请对该族用更细的网格复核
   （先 `python ../E2_transmutation_atlas/e2a_atlas.py --hq --h_scale 0.75 --tag hq_fine --samples <样品> --ms <m>`，
   再 `python e3_stability.py --hq --h_scale 0.75 --tag hq_fine --samples <样品> --ms <m> --zmax 2000`），并在 Progress.md 中报告。

## 上传清单
`nonlocal_modes/results/E3_stability/main/` 或 `hq/`（以及 `hq_coarse`）整个目录。Progress.md 中记录：每个 (样品, m) 稳定的 b 区间；任何 m≥3 的稳定段；传播与线性谱不一致的情形。

---

# E3b —— 低功率失稳窗口：阈值、失稳模式的结构、按角向通道分解的谱（第二轮，2026-10）

**为什么要做**（Progress.md C8、开放问题 Q6）：hq 全扫描发现 m=2 在**所有**样品（含圆盘）中都有低功率振荡失稳窗口
（圆盘 b_c≈0.53），三角形把窗口拉宽到 b_c≈2.7，三角形中 m=1 也失稳。E3 的 b 间隔 0.2、只记录增长率，
回答不了“窗口边界在哪里、失稳的是什么模式、为什么三角形特别宽”。E3b 加密 b 网格、在交界处二分、对每个不稳定点做角向谐波分析。

**机理假说**（待本实验检验）：Bogoliubov 微扰 δq = u e^{λξ} + v* e^{λ*ξ}，电荷 m 的涡旋中 u ∝ e^{i(m+L)θ}、v ∝ e^{i(L−m)θ}，
L 为微扰的角向指数。圆盘中不同 L 互不耦合；C_n 边界把 L 与 L±n 耦合，两支 Krein 符号相反的频率相遇处产生新的失稳窗口。
κ = L mod n（取在 (−n/2, n/2]）：κ=0 的微扰保持 C_n 对称，κ≠0 的微扰破坏 C_n 对称（E4 中 μ 权重漂移的来源）。

## 运行
```bash
cd nonlocal_modes/experiments/E3_stability
python e3b_instability_windows.py --quick                          # 约 10 分钟 (圆盘与三角形 m=2, h=0.4)
python e3b_instability_windows.py --hq --cases poly3:2 --job poly3_m2       # 一个情形
python e3b_instability_windows.py --hq --cases disk:2 --channels --job disk_m2_channels   # 圆盘通道分解谱
# 网格收敛表与旋转检验 (同一组 b 点上的细网格, 不做传播):
python e3b_instability_windows.py --hq --h_scale 0.5 --tag hq_h005 --cases poly3:2 --b_range 2.2 5.0 --spec_step 0.4 --no_evolution --job poly3_m2_refine
python e3b_instability_windows.py --hq --cases poly3:2 --b_range 1.8 3.4 --theta0 15 --no_evolution --job poly3_m2_rot15_refine
python e3b_instability_windows.py --hq --table_only                # 拆分运行结束后重建 e3b_summary.csv
```
超算上不需要手打这些命令：`python hpc/make_jobs.py --phase 2` 会生成全部第二轮任务（见 [`../hpc/README.md`](../hpc/README.md) 第 7 节）。

情形（`--cases`，写法 `样品:m`）：rect:2、poly3:2、poly4:2、poly5:2、poly6:2、disk:2、poly3:1、disk:1。
b 范围默认 0.05–1.4（三角形 m=2 到 5.0，正方形 m=2 到 1.6），hq 档 b 步长 0.05、每 0.1 算一次谱（`--spec_step` 可改），交界处二分 4 次（阈值精度约 0.006）。

## 输出与预期结果
* `e3b_<样品>_m<m>_spectrum.csv`：每个 b（含二分加密的点）的增长率、最不稳定本征值、`L`、`kappa`、u/v 的前三个谐波（`top_u`、`top_v`）、`v_over_u`；
* `e3b_<样品>_m<m>_growth.png`：增长率 vs b，每个不稳定点标出 L；
* `e3b_<样品>_m<m>_mode_b<值>.png`：每个窗口最不稳定点与最靠近阈值的点的模式图（|u|、arg u、|v|、扰动后的光强）；
* `e3b_<样品>_m<m>_evolution_b<值>_*.png/.csv`：沿最不稳定模加 1e-3 扰动后的传播（非线性结局：分裂、振荡还是饱和）；
* `e3b_disk_m<m>_channels.csv/png`（`--channels`）：圆盘中 Im λ ≤ 1 的全部本征值按 L 分类（颜色）并标出 Krein 符号（实心圆 +、空心三角 −）；
* `e3b_summary.csv`：每个情形一行（阈值区间 `thresholds`、窗口 `windows`、各窗口最大增长率与 L/κ）。

**预期**（云端 h=0.25–0.4 探索）：
* 圆盘 m=2：窗口 b < 0.5，L = −2（u ∝ e^{0iθ}、v ∝ e^{−4iθ}）——双荷涡旋分裂为两个单位涡旋的四极模（与 BEC 双量子化涡旋的分裂失稳同类）；
* 三角形 m=2：b≈1.5 处为 L = −1 与 +2 混合（相差 n=3，κ = −1）；低 b 处为 L = −2 与 +1 混合（κ = +1）；
* 三角形 m=1：L = −2 与 +1 混合（κ = +1）；
* 五边形、六边形 m=2：窗口与圆盘相同（只有高阶耦合）；
* 通道谱：三角形中出现窗口的 b 处，圆盘的 L 与 L±3 两支（Krein 符号相反）应该相交或很接近。

**已知问题（请特别关注）**：三角形 m=2 的窗口上端随网格明显移动（h=0.4: 1.9–2.0；h=0.15: 2.4–2.6；h=0.1: 2.6–2.8）。
云端定点诊断（2026-10-12）：窗口上部（b ≳ 2.2）是 **κ=0 的弱振荡失稳**（L=−3 与 L=0 混合，u 的主导谐波 l=−1、+2，即中心 −1 与外层 +2
两个分量之间的“卫星环呼吸”模，Im λ≈0.66–0.76），它的增长率与 h=0.1 的离散误差同量级：b=2.8 处 h=0.1/0.075/0.05 的增长率为
稳定/2.1e-3/4.3e-3（Im λ 都是 0.764），由后两点按 g = g₀ − c h² 外推 g₀ ≈ 6e-3（连续极限下 b=2.8 很可能仍不稳定；
外推式在 h=0.1 处为负，即本征值回到虚轴，与 h=0.1 判稳定一致）；b=3.6 处三档网格都稳定，所以连续极限的上端在 2.8–3.6 之间。
所以“阈值区间”在这里没有意义，改为**在同一组 b 点（2.2, 2.6, …, 5.0）上算四档网格**（h=0.15/0.1/0.075/0.05，
tag `hq_coarse`/`hq`/`hq_h0075`/`hq_h005`），逐点按 h² 外推增长率，得到连续极限下的失稳上端；另做多边形相对网格转 15° 的旋转检验
（文件名带 `_rot15`）：若 κ=0 模的增长率随转角明显改变，说明网格 C_4 各向异性也在起作用。

**耗时（hq 单核，云端实测每个 shift：h=0.1 13 s、h=0.075 30 s、h=0.05 100 s）**：每个谱点 29 个 shift，h=0.1 约 6 分钟；
三角形 m=2（b 到 5.0）约 7–8 h，其余情形 2–3 h；h=0.05 加密约 11 h（内存峰值约 5 GB），h=0.075 约 4 h；通道谱约 1 h。

**上传**：`results/E3_stability/hq/` 中新增的 `e3b_*` 文件与 `log_*.txt`，以及 `hq_coarse/`、`hq_h0075/`、`hq_h005/` 中的 `e3b_*` 与日志。
Progress.md 中记录：每个情形的阈值区间与窗口内的 L/κ；三档网格 + 旋转检验的三角形阈值；通道谱中 L 与 L±3 两支是否在窗口处相交。
