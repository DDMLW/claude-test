# E1 —— 对称性传递定律与保角半径普适性（论文 Fig.1）

**物理问题**：材料完全均匀，只有截面是正 n 边形。边界的 C_n 对称性以什么形式进入光束感受到的非线性折射率？

**理论预言**（docs/theory.md 第 2 节）：
1. 中心处 C_n 对称光束（总能流 U）产生的折射率中，最低阶各向异性项为 δn_n(r) cos nθ，
   **δn_n(r) = U/(π n(n+1)) · (r/ρ)^n**，ρ 为保角半径；与光束形状无关，只出现 k = 0, n, 2n… 次角谐波；
2. 保角半径相同的所有正多边形与圆盘，C_n 对称态的 U(b) 曲线几乎完全重合；
3. m=1 涡旋（n≥3）受对称性保护：中心荷始终为 1，但亮环被 C_n 调制，调制深度随 n 急剧减小。

## 运行
```bash
cd nonlocal_modes/experiments/E1_symmetry_transfer
python e1a_index_anisotropy.py --quick     # 约 1 分钟
python e1b_protected_modes.py --quick      # 约 1 分钟
python e1a_index_anisotropy.py             # 正式, 约 10 分钟
python e1b_protected_modes.py              # 正式, 约 30-60 分钟
```

## 输出与预期结果
### e1a_index_anisotropy.py
* `e1a_power_law.png`：四个子图 (n=3,4,5,6)，点 = 数值，线 = 理论，log-log 斜率 = n。**点应落在线上。**
* `e1a_index_anisotropy.csv`：
  * `part=power_law`：`ratio` 列（数值/理论）应为 1.0000 ± 1e-3（n=4 在 r=1 处可能偏离几个百分点：那里信号极小，网格误差占比大，属正常）；
    `max_forbidden`（非 n 倍数的谐波）应比 `dn_n` 小 3 个量级以上；
  * `part=universality`：不同宽度的高斯光束、环形光束，`ratio` 都应 = 1.0000；
  * `part=linearity`：U 从 2 到 40，`ratio` 都应 = 1.0000；
  * `part=rectangle`：矩形的数值系数（正方形应为 0）。

### e1b_protected_modes.py
* `e1b_families.png`：左图 U(b)——所有样品 (poly3–poly6, disk, ρ=8) 的 m=0 曲线重合、m=1 曲线重合；右图 m=1 环上 cos nθ 调制深度 vs U（log-log），n 越大越小。
* `e1b_<样品>_m<m>_summary.csv`：每个解的诊断量。检查 m=1 的 `central_charge` 都 = 1，`mu_weight` ≈ 1。
* `e1b_<样品>_m<m>_b2.00.png`：b=2 的模式图，m=1 涡旋在三角形中亮环可见三重调制。
* 预实验数值（h=0.4）：b=0.5 时 m=1 调制深度 n=3: 1.5e-2, n=4: 1e-3, n=5: 1.3e-4, n=6: 1.9e-5。

## 上传清单
`nonlocal_modes/results/E1_symmetry_transfer/main/` 整个目录。Progress.md 中记录：e1a 的 ratio 范围；e1b 中不同样品 U(b) 的最大相对差。
