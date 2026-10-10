# E5 —— 热边界可重构调控（论文 Fig.3 / Fig.6）

**思想**：不改样品形状，而是让圆形样品的**边界温度**按 Θ cos n(θ−θ₀) 分布（实验上可用沿圆周排布的加热片/珀尔帖元件）。
热传导方程是线性的，样品内多出一个与光功率无关的调和折射率 Θ(r/R)^n cos n(θ−θ₀)：
各向异性强度 Θ 可连续调节、可反号、可旋转，还可以沿传播方向分段改变。

## E5a —— 稳态：r_s(Θ) 与取向翻转
```bash
cd nonlocal_modes/experiments/E5_thermal_reconfiguration
python e5a_boundary_heating.py --quick     # 约 5 分钟 ((3,2), (4,3) 两种情形, Θ≤1)
python e5a_boundary_heating.py             # 正式, 约 30-60 分钟, 六种 (n,m)
```
输出：`heating_n<n>_m<m>.csv`（每个 Θ：中心荷、卫星数、r_s 检测值/拟合值、第一颗卫星方位角）、`heating_rs.png`（r_s vs Θ 的 log-log 图与拟合斜率）、`heating_mode_*.png`。

预期（docs/theory.md 第 3.5 节）：
| (n, m) | 中心荷 | 卫星 | r_s ∝ Θ^s |
|---|---|---|---|
| (3, 2) | −1 | 3 | s = 1 |
| (4, 3) | −1 | 4 | s = 1/2 |
| (3, 3) | 0 | 3 | s = 1/3 |
| (5, 3) | −2 | 5 | s = 1 |
| (6, 4) | −2 | 6 | s = 1/2 |
| (2, 2) | 0 | 2 | s = 1/2（四极加热 → 可调“双涡旋”） |

* Θ 反号 → 卫星环整体转 π/n（看 `sat_angle` 列）；
* 预实验（h=0.4/0.2）：(3,2) 拟合斜率 **1.004**；(4,3) 0.446（理论 0.5），Θ 反号卫星环转 45°；(6,4) 中心 −2、6 个卫星，r_s ∝ Θ^0.52；(5,3) 中心 −2 在 Θ ≥ 2 时出现。
* n=5、6 的嬗变规则主要靠本实验展示（多边形样品中太弱，见 E0d）。

## E5b —— 动态控制：关断、切换、旋转
```bash
python e5b_dynamic_control.py --quick      # 约 5 分钟
python e5b_dynamic_control.py              # 正式, 约 1-2 小时 (xi 到 300)
python e5b_dynamic_control.py --Theta 2 --tag T2    # 更强的加热
```
起点：C_3 加热下 m=2 的嬗变定态（中心 −1 + 3 卫星）。四种情形：
* `switch_off`：ξ > ξ_s 关断加热 → 预期卫星收回、中心荷回到 +2（**可逆嬗变**；预实验已观察到）；
* `switch_C3_to_C4`：换成 C_4 图样 → m=2 在 C_4 中是 μ=n/2 简并情形，记录结构如何演化；
* `rotate_slow` / `rotate_fast`：图样以 Ω 旋转 → 慢转时卫星应绝热跟随（`tracks.png` 中方位角线性增长），快转时跟不上。

输出：`control_start_state.png`、`control_<case>_snapshots.png`、`control_<case>_tracks.csv/png`（奇点半径/方位角随 ξ）、`control_summary.csv`。

注意：横向热传导按稳态处理，要求图样沿 ξ 的变化远慢于横向热扩散；E5b 的旋转部分主要是拓扑动力学的理论展示，实验可行性需在论文中讨论。

## 上传清单
`nonlocal_modes/results/E5_thermal_reconfiguration/main/`（及 `T2` 等）。Progress.md 中记录：各 (n,m) 的拟合斜率；Θ 反号时的方位角变化；E5b 四种情形的末态中心荷与卫星轨迹的定性描述。
