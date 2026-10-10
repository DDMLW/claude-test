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
```
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

**判读规则**：
1. 线性谱判不稳定但传播中 `survived=True` 且偏离度有界振荡：记为“弱不稳定、非线性饱和”（不是数值错误，见 Progress.md 中表面偶极的例子）；
2. 线性谱判稳定但传播破坏：极少见，可能是 shift 扫描漏掉了本征值，请用 `--omega_max 5` 重算该族并报告；
3. 若发现某个 m ≥ 3 的嬗变态在一段 b 内**稳定**，这是重要结果：请对该点用 `--h 0.15 --tag h015`（需先用同样参数跑 E2a 该样品）复核，并把传播距离加到 `--zmax 2000`。

## 上传清单
`nonlocal_modes/results/E3_stability/main/` 整个目录。Progress.md 中记录：每个 (样品, m) 稳定的 b 区间；任何 m≥3 的稳定段；传播与线性谱不一致的情形。
