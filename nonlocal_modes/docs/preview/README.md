# 参考预览图（不是正式结果）

这些图是 Claude 在云端用 **quick 模式（粗网格 h = 0.3–0.4、少量参数点、短传播距离）** 跑出来的，
只用来让你知道“正确运行时图大概长什么样”。正式结果以你在本地用默认参数（h = 0.2）跑出的 `results/<实验>/main/` 为准，
细节（奇点位置、阈值、增长率）会因网格加密而改变；粗网格下卫星环小于 2h 的结构本来就分辨不出来。

| 文件 | 对应实验 | 看什么 |
|---|---|---|
| `E0a_convergence.png` | E0a | 左：Poisson 误差二阶下降；右：双涡旋 D_s 随 h 收敛 |
| `E0b_twin_vortex_m2_b2.png`, `..._m3_b4.64.png` | E0b | 导师 2010 年 Fig.1：奇点沿短轴 x 排列 |
| `E0b_fig2b_fig3a.png` | E0b | W_η ≥ W_ζ；D_s 随 d_η 的变化 |
| `E0c_surface_vortex_b3.png` | E0c | 导师 2008 年 Fig.2：贴着绝热界面的非正则涡旋 |
| `E1a_power_law.png` | E1a | 点（数值）落在线（理论 U/(πn(n+1))(r/ρ)^n）上 |
| `E1b_families.png` | E1b | 左：不同多边形（ρ 相同）的 U(b) 重合；右：m=1 涡旋的 C_n 调制随 n 急剧减小 |
| `E2a_atlas_phase.png`, `E2a_atlas_amplitude.png` | E2a | 嬗变周期表：红圈 +1、青圈 −1，中心数字为中心荷 |
| `E2b_scaling.png` | E2b | r_s vs ρ 的 log-log 斜率与理论比较 |
| `E3_poly3_m2_stability.png`, `E3_poly3_m3_stability.png` | E3 | 三角形中 m=2 嬗变态稳定、m=3 不稳定 |
| `E4_poly3_m2_worldlines.png`, `E4_poly3_m2_snapshots.png` | E4 | 正则 m=2 光束射入三角形样品：中心荷（黑线）由 +2 变为 −1 |
| `E5a_heating_rs.png` | E5a | r_s vs 加热幅度 Θ 的 log-log 斜率 |
| `E5b_switch_off_tracks.png` | E5b | 关断加热后卫星收回中心（嬗变可逆） |
