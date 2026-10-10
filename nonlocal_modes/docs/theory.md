# 理论推导

本文档给出研究方案中各条理论预言的推导。代码 `nlmodes/theory.py` 实现了其中可直接算出数字的公式，
实验脚本会把数值结果与之逐项对比。**凡是标注“预言”的结论都必须以数值结果为准**；若数值与理论不符，
先查数值（网格、分辨率），再修正理论，并记录在 Progress.md 中。

---

## 1. 模型（与导师论文完全相同）

$$
i\,\partial_\xi q=-\tfrac12\Delta_\perp q-(n_{\rm lin}+n)\,q,\qquad
\nabla\cdot(\kappa\nabla n)=-|q|^2\ \ (\text{热介质内}),
$$

边界：$n=n_b$（热沉，常取 0；E5 中为边界温度图样）或 $\partial_\nu n=0$（绝热界面），光场在样品外为 0。
定态 $q=w\,e^{ib\xi}$；能流 $U=\iint|q|^2$。

**物理单位**（同导师 PRA 2008）：$q=(k_0^2r_0^4\alpha\beta/\kappa n_0)^{1/2}A$，$n=k_0^2r_0^2\,\delta n/n_0$，
横向坐标以 $r_0$、纵向以 $L_{\rm dif}=k_0r_0^2$ 为单位，于是 $U=(k_0^2r_0^2\alpha\beta/\kappa n_0)\,P$。
例：铅玻璃 $n_0=1.8,\ \beta\sim10^{-5}\,{\rm K^{-1}},\ \alpha\sim0.01\,{\rm cm^{-1}},\ \kappa\sim1\,{\rm W\,m^{-1}K^{-1}}$，
$\lambda=0.5\,\mu$m。**注意**：$r_0=50\,\mu$m 时 $L_{\rm dif}\approx3$–6 cm（取决于 $k_0$ 是否含 $n_0$），
厘米级样品只对应 $\xi\sim1$–3；$r_0=20\,\mu$m 时 $L_{\rm dif}\approx0.5$–1 cm，对应 $\xi\sim5$–10。
因此 E4 中“嬗变发生的传播距离”是判断实验可行性的关键量。

---

## 2. 对称性传递定律

### 2.1 保角映射表示的格林函数
设热介质占据单连通区域 $\Omega$（边界恒温），$f:\Omega\to$ 单位圆盘为保角映射，$f(0)=0$，$f'(0)=1/\rho>0$
（$\rho$ 称为 $\Omega$ 关于原点的**保角半径**）。$-\Delta G=\delta$ 的 Dirichlet 格林函数为

$$
G(z,z')=\frac1{2\pi}\ln\left|\frac{1-f(z)\overline{f(z')}}{f(z)-f(z')}\right|,\qquad n(z)=\int_\Omega G(z,z')\,|q(z')|^2\,d^2z'.
$$

### 2.2 正 n 边形的 Schwarz–Christoffel 映射
顶点指向 $\theta=0$ 的正 $n$ 边形（$n\ge3$）：
$$z=F(w)=\rho\int_0^w(1-s^n)^{-2/n}ds=\rho\Big[w+\tfrac{2}{n(n+1)}w^{n+1}+O(w^{2n+1})\Big],$$
外接圆半径 $R_c=\rho\,I_n$，$I_n=\int_0^1(1-s^n)^{-2/n}ds=\dfrac{\Gamma(1/n)\Gamma(1-2/n)}{n\,\Gamma(1-1/n)}$
（$I_3=1.7666,\ I_4=1.3110,\ I_5=1.1745,\ I_6=1.1129$；代码 `geometry.polygon_integral_In`）。反函数
$$f(z)=\frac z\rho\Big[1+c_n\Big(\frac z\rho\Big)^n+O\big((z/\rho)^{2n}\big)\Big],\qquad c_n=-\frac{2}{n(n+1)}.$$

### 2.3 展开与定律
对 $|z|,|z'|\ll\rho$：
$$\ln|f(z)-f(z')|=\ln|z-z'|-\ln\rho+\mathrm{Re}\Big[\frac{c_n}{\rho^n}\sum_{j=0}^{n}z^jz'^{\,n-j}\Big]+\dots,\qquad
\ln|1-f(z)\overline{f(z')}|=-\mathrm{Re}\frac{z\bar z'}{\rho^2}+\dots$$
记光强的复矩 $M_k=\int z'^k|q|^2d^2z'$。若光强关于中心 C_n 对称，则 $M_k=0$（除非 $n\mid k$），于是
$$
n(z)=n_{\rm free}(z)+\frac U{2\pi}\ln\rho+\underbrace{\frac{U}{\pi n(n+1)}\Big(\frac r\rho\Big)^n\cos n\theta}_{\text{最低阶各向异性}}+\dots
$$
其中 $n_{\rm free}=-\frac1{2\pi}\int\ln|z-z'|\,|q|^2$ 是无界介质的响应。由此得到两个结论：

1. **对称性传递定律**：最低阶各向异性系数 $U/(\pi n(n+1))\cdot\rho^{-n}$ 只依赖总能流与保角半径，与光束形状无关；
   角谐波只出现 $k=0,n,2n,\dots$。（E1a 验证：数值/理论 = 1.0000。）
2. **保角半径普适性**：各向同性部分 = 自由响应 + 常数 $(U/2\pi)\ln\rho$；后者只平移传播常数。
   因此保角半径相同的任何正多边形（及半径为 $\rho$ 的圆盘）中，C_n 对称态的 $U(b)$ 曲线在 $O((r/\rho)^n)$ 精度内重合。（E1b 验证。）

矩形（C_2）没有简单闭式，E1a 给出数值系数。

---

## 3. 角赝动量与选择定则

### 3.1 守恒律
设样品（以及边界温度图样）在绕中心转 $2\pi/n$ 的操作 $\hat R$ 下不变。若场满足 $\hat R q=e^{-2\pi i\mu/n}q$
（不可约表示 $\mu\in(-n/2,n/2]$），则 $|q|^2$ C_n 对称 → 折射率 C_n 对称 → 方程 C_n 等变，
所以 **μ 在传播中严格守恒**，定态也可按 μ 分类。电荷 $m$ 的正则涡旋 $\propto e^{im\theta}$ 属于 $\mu=m\bmod n$。
（E4 验证：μ 权重的漂移只来自网格对 C_n 的微小破坏。）

### 3.2 谐波结构
属于 μ 的场只含角谐波 $l\equiv m\ (\mathrm{mod}\ n)$：$w=\sum_j w_{l_j}(r)e^{il_j\theta}$，$l_j=m-jn$（及 $m+jn$），
在中心正则性要求 $w_l(r)\sim a_l r^{|l|}$。把各向异性 $\varepsilon\,r^n\cos n\theta$ 当作对正则涡旋的微扰，
每作用一次角动量改变 $\pm n$（线性化同时包含 $w$ 与 $\bar w$，不改变这个计数），故 $a_{l_j}=O(\varepsilon^j)$。

### 3.3 中心荷与“牛顿多边形”
在半径 $r$ 处 $|a_{l_j}|r^{|l_j|}\sim\varepsilon^j r^{|l_j|}$ 最大的谐波主导相位。令 $t=\ln(1/r)$、$e=\ln(1/\varepsilon)$，
主导项使 $je+|l_j|t$ 最小。随 $r\to0$，主导谐波沿点集 $\{(|l_j|,j)\}$ 的**下凸包**从 $(|m|,0)$ 移向 $|l|$ 最小的点：

* **中心荷** $v$ = 凸包终点的 $l$（$|l|$ 最小、同 $|l|$ 时取 $j$ 最小者）；
* 凸包每条边 $(|l_a|,j_a)\to(|l_b|,j_b)$ 对应一个卫星环：环上奇点总荷 $l_a-l_b$（单位荷，个数为 $n$ 的整数倍，C_n 对称排列），
  环半径 $r\propto\varepsilon^{s}$，$s=(j_b-j_a)/(|l_a|-|l_b|)$。

例（代码 `theory.selection_rule`）：

| (n, m) | μ | 中心荷 v | 卫星 | 环半径标度 |
|---|---|---|---|---|
| (2, 2) 矩形 | 0 | 0 | 2 个 +1 | $\varepsilon^{1/2}$（导师的双涡旋） |
| (2, 3) 矩形 | 1 | +1 | 2 个 +1 | $\varepsilon^{1/2}$ |
| (3, 2) | −1 | −1 | 3 个 +1 | $\varepsilon^{1}$ |
| (3, 3) | 0 | 0 | 3 个 +1 | $\varepsilon^{1/3}$ |
| (4, 3) | −1 | −1 | 4 个 +1 | $\varepsilon^{1/2}$ |
| (4, 4) | 0 | 0 | 4 个 +1 | $\varepsilon^{1/4}$ |
| (5, 3) | −2 | −2 | 5 个 +1 | $\varepsilon^{1}$ |
| (6, 4) | −2 | −2 | 6 个 +1 | $\varepsilon^{1/2}$ |
| (3, 5) | −1 | −1 | 3+3 个 +1（两个环） | $\varepsilon^{1/3}$, $\varepsilon^{1}$ |
| (n, m<n/2) | m | m | 无（受保护，但非正则） | — |

**超敏感性**：$s=1/(|m|-|v|)$ 可以很小（如 (4,4) 的 1/4），高阶涡旋对极弱的各向异性也会产生宏观分裂——
这也意味着名义上圆形、实际有微小缺陷的样品中高阶涡旋同样会分裂（E0d 中网格自身的 C_4 误差正是这样一个“缺陷”）。

### 3.4 简并情形与存在阈值
若 $n\mid 2m$，谐波 $-m$ 出现在第 $j=2m/n$ 阶，而 $\pm m$ 在各向同性极限下简并 → 共振耦合，简单微扰论失效：
线性极限下的本征态是实驻波（$\cos m\theta$、$\sin m\theta$ 型，分属 C_nv 的不同一维表示），
复涡旋态需要非线性锁相，**预期只在功率（或 b）阈值以上存在**，阈值随各向异性增强而升高。
矩形中 m=2 双涡旋在大纵横比时存在下阈值（导师 2010 年已观察到）即属此类。正方形中 m=2（μ=n/2）中心荷为 +2。

### 3.5 各向异性强度 ε 的两种来源
* 形状：$\varepsilon_{\rm shape}=U/[\pi n(n+1)\rho^n]$（$r^n\cos n\theta$ 的系数，第 2 节）；
* 边界加热（E5）：边界温度 $\Theta\cos n(\theta-\theta_0)$ 在圆盘内产生调和场 $h=\Theta(r/R)^n\cos n(\theta-\theta_0)$，
  $\varepsilon_{\rm heat}=\Theta/R^n$，与光功率无关，可正可负、可旋转。等效温度 $\Theta_{\rm eq}=UR^n/[\pi n(n+1)\rho^n]$。

于是**固定 U 时 $r_s\propto\rho^{-ns}$**（E2b），**固定 U 时 $r_s\propto\Theta^{s}$**（E5a），Θ 反号 ⇔ 图样转 $\pi/n$ ⇔ 卫星环转 $\pi/n$。

**极值原理**：调和函数在区域内部没有极大值，外加的 $h$ 本身不能导引光；只有与自陷的热透镜叠加才形成定态。
对 n=2（四极加热）中心是鞍点，Θ 过大时光束会被推离中心（E5a 中 n=2 的 Θ 上限取小）。

---

## 4. 稳定性（定性预期）
* 线性化算子是哈密顿型，本征值 $(\lambda,-\lambda,\bar\lambda,-\bar\lambda)$ 成组；离散化保持对称性，因此不会产生虚假增长率。
* 已知：圆柱样品中 $m\le2$ 稳定（Kartashov 2007），矩形中 m=2 稳定、m≥3 不稳定（Ye 2010）。
* 嬗变态的稳定性是开放问题：C_n 各向异性会使原本简并的方位扰动通道失谐，可能抑制（也可能加剧）不稳定性。
  预实验：三角形中 m=2 嬗变态稳定，m=3 振荡不稳定。完整结论见 E3。
* 线性谱判为弱不稳定、但传播中扰动非线性饱和（长期有界振荡）的情形要单独标注（E0c 中表面偶极即属此类，见 Progress.md）。

---

## 5. 旋转边界图样（E5b）
在随图样转动的坐标系 $\theta'=\theta-\Omega\xi$ 中，$i\partial_\xi Q=-\tfrac12\Delta Q-V(r,\theta')Q-\Omega L_zQ$，$L_z=-i\partial_\theta$。
小 Ω 时卫星绝热跟随（世界线为螺旋），Ω 超过内部模式频率差时跟不上（动力学对称性恢复）。
注意：横向热传导按稳态处理，要求图样沿 ξ 的变化尺度远大于样品横向尺寸（物理单位下的适用性需讨论）。
