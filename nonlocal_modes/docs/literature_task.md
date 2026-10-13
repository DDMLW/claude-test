# 投稿前文献调研任务书（第六轮，交给师弟；2026-10-13 布置）

> 本任务**全部由师弟在本地完成**（云端环境连不上 arXiv / OpenAlex / Semantic Scholar，见 references.md F 节说明）。
> 不改代码、不改结果；只更新文献文件，并在 `Progress.md` §6 记一行运行记录。
> 预计 5–7 个工作日。发现“撞车”（别人已经做了同一件事）时**立即**告诉师兄，附 PDF，不要等全部做完。

## 0. 我们要投的论文到底声称了什么（逐条检验的对象）

| 编号 | 声称（新意） | 本文依据（Progress.md） |
|---|---|---|
| N1 | 材料完全均匀，离散旋转对称性只来自**样品边界**；热非局域响应把边界的 C_n 对称性按解析定律 δn = U/(πn(n+1))·(r/ρ)^n cos nθ 传到光束处；多边形与等保角半径 ρ 的圆盘 U(b) 重合 | C1、C2 |
| N2 | 由边界诱导的涡旋**拓扑荷嬗变选择定则**（中心荷 = m − jn 中绝对值最小者、卫星数、角赝动量 μ）与卫星半径的**牛顿多边形分数标度律** r_s ∝ ε^s | C3–C5 |
| N3 | m=2 涡旋的低功率失稳窗口被边界对称性拉宽（三角形把 b_c 从 0.50 推到 ≈2.97）；机理 = 圆盘 Bogoliubov 通道的 **Krein 碰撞 + C_n 选择定则 κ = L mod n**，Hamilton–Hopf 判据 \|Δω\| = 2g 定量预言窗口端点 | C8、C12、C13 |
| N4 | 动态嬗变的可观测窗口：LG 光束在 ξ≈2–12 完成嬗变，稳定区保持到 ξ=600 | C9 |
| N5 | **四极边界加热的可调“双涡旋”**：m=2 分裂为两个 +1 涡旋，半间距 r_s ∝ Θ^{1/2}；Θ_c = 1.176 处破坏 C_2 的叉式分岔（增长率三位数字预言） | C10 |
| N6 | **受迫分裂模图像**：关断后以 L=−2 分裂模频率的一半自由共转；旋转热图样拖动双涡旋时间距按 ω_rot/(ω_rot − 2Ω) 放大（无可调参数）；同向锁定受一阶 / 二阶共振限制，反向锁定受破坏 C_2 的平移模在旋转坐标系中频率过零（\|Ω\| ≈ ω_odd）限制 | C10（第五轮 E5e） |

## 1. 检索词（中英文都要查；时间不限，但 2020–2026 年要逐篇看摘要）

* **A 边界 / 样品形状**：`thermal nonlocal` AND (`boundary` OR `sample shape` OR `cross-section` OR `polygon` OR `square` OR `rectangular` OR `triangular`) AND (`vortex` OR `soliton`)；`boundary force soliton nonlocal`；`lead glass` `vortex` `boundary`；中文：热非局域 / 强非局域 + 边界 / 样品形状 / 截面 + 涡旋孤子
* **B 涡旋嬗变 / 离散对称**：`vortex transmutation`；`angular pseudomomentum`；`discrete rotational symmetry` `vortex soliton`；`topological charge` `selection rule` `discrete symmetry`；`off-axis singularities` `charge transmutation`
* **C 双量子化涡旋分裂 / Krein**：`doubly quantized vortex` `splitting` `instability`；`doubly charged vortex` `anisotropic trap` / `elliptic trap` / `quadrupole deformation` `splitting`；`Krein signature` `vortex`；`Hamiltonian Hopf` `vortex soliton`；`negative energy mode` `vortex`
* **D 旋转图样 / 旋转晶格 / 旋转势**：`rotating lattice` `soliton` `critical angular velocity`；`rotary soliton` `Bessel lattice`；`rotating potential` `vortex dragging`；`rotating anisotropic trap` `vortex pair`；`azimuthon`；`rotating dipole soliton nonlocal`
* **E 热边界的主动调控**：`boundary temperature` / `heating pattern` / `thermo-optic control` + `soliton steering` / `reconfigurable`；`optically induced thermal` `nonlocal` `control`
* **F 标度律**：`vortex splitting` `scaling` `perturbation` `anisotropy`；`Newton polygon` `vortex`（大概率没有，查一下即可）

## 2. 必须逐条回答的问题（按优先级）

1. **KVT 2007（references.md B9，Opt. Express 15, 9378）**：他们报告圆柱热介质中 m ≤ 2 稳定、“与样品半径无关”。我们的圆盘 m=2 在**低功率端**有窄失稳窗口
   （L=−2 分裂模，b_c R² = 31.8；ρ=8 时 b < 0.50，R=10 时 b < 0.32；h=0.05 收敛），之上稳定。请读他们的稳定性图：b（或能流 U）的扫描下限是多少？
   是否覆盖 b R² < 32？若覆盖且判为稳定，请找出模型差别（方程、边界条件、归一化）。他们的“与半径无关”就是我们的严格标度律 b → bR²。**这一条最重要。**
2. **BEC 中的对应物**（C 组检索）：有没有人在**各向异性（椭圆 / 四极形变）势阱**中让双量子化涡旋分裂成两个间距 ∝ ε^{1/2} 的单涡旋，
   并用**旋转的各向异性势阱**拖动这对涡旋、讨论锁定的临界角速度？这是 Fig.6（N5、N6）最可能的先行工作，务必仔细查
   （线索：Möttönen 2003、Kawaguchi & Ohmi 2004、Huhtamäki 2006、Kuopanportti、Lundh & Nilsen 2006、Kapitula–Kevrekidis 2013）。
3. **导师组的双涡旋**（Ye 等 OL 35, 628 (2010)，A2）：他们的双涡旋是怎样定义的（矩形样品中的 m=2？还是两个涡旋的束缚态？），我们的“C_2 边界 / 四极加热诱导的双涡旋”
   如何引用、哪些是新的（Θ 可调、Θ^{1/2}、叉式分岔、旋转拖动）。
4. **多边形 / 异形截面热介质 2020–2026**（A 组）：references.md B12（七边形热介质中孤子轨道进动，Int. J. Theor. Phys. 2025）的作者与同组其他工作；
   国内课题组（如华南师大郭旗组、浙江师大 / 广东相关组等强非局域方向）关于矩形 / 椭圆 / 多边形截面、边界条件对孤子影响的论文——是否有人研究过涡旋的拓扑荷？
5. **旋转晶格 / 旋转势中的锁定**（D 组）：Kartashov–Vysloukh–Torner 2004 “Rotary solitons in Bessel optical lattices”（PRL 93, 093904，待核对）、
   Sakaguchi–Malomed 2007、Cuevas 等 2007（F31–32）：他们给出的临界角速度由什么机理决定？有没有用内部模的共振（ω/2、ω/4）或 C_2 破缺模频率过零解释？
6. **Ferrando 等的涡旋嬗变**（C15–16、C36）：他们的选择定则（m − kn）与我们的中心荷规则是否完全一致；有没有人做过“由样品边界（而不是晶格）诱导的嬗变”。
7. **旋转孤子 / azimuthon 与本文的区别**（D 组、F33–34）：我们的双涡旋在旋转图样中是**被迫同步**的定态（旋转坐标系中的定态），不是自发旋转的 azimuthon——
   请确认文献中有没有“外加旋转边界 / 旋转加热”驱动的非局域孤子。
8. **热边界主动调控**（E 组）：有没有人用分段 / 图样化的边界加热控制涡旋或孤子簇（不仅是基模的偏转）。
9. **分数标度律**（F 组）：有没有人推导过卫星涡旋位置随各向异性强度的分数幂律。
10. **核对 references.md A–F 全部“待核对”条目**：作者、标题、期刊、卷、页、年份、DOI。

## 3. 交付物

1. **`nonlocal_modes/docs/references.md`**：每条核对后在行首加 ✓（或改正后加 ✓）；新增的相关文献放到新的 **G 节**（每条附 DOI 与一句话说明与本文的关系）。
2. **`nonlocal_modes/docs/literature_review.md`**（新建）：
   * 一张表：N1–N6 每一条 → 最接近的 2–4 篇先行工作 / 他们做了什么 / 与我们的区别 / 风险（无 / 部分重叠 / 撞车）；
   * 第 2 节 10 个问题的逐条回答（每条一段，带引用）；
   * 一段新意陈述（中文即可，师兄改成英文放进引言）；
   * 期刊建议：Chaos, Solitons & Fractals；Communications in Nonlinear Science and Numerical Simulation；Physical Review A；Photonics Research；
     Optics Express；Physica D——各自的**中科院分区（最新版）**、影响因子、近两年是否发表过同类工作（列 1–2 篇）。目标是中科院一区或以上。
3. **`nonlocal_modes/paper/refs.bib`**：核对过的 BibTeX（键名用 `第一作者年份关键词`，如 `Ye2010twin`、`Kartashov2007cyl`），供正文草稿引用。
4. 在 `Progress.md` §6 记一行（日期、做了什么、找到几篇相关 / 有无撞车），§4-13（KVT 2007）写上结论。

## 4. 方法建议

* Google Scholar、Web of Science、Scopus、arXiv、CNKI（中文期刊与学位论文也看，例如《物理学报》《光学学报》）。
* 本仓库自带的 `paper-search-pro` / `research-lit` skill（见仓库根目录 README：在能连 OpenAlex / Semantic Scholar / arXiv 的 Claude Code 会话里用，
  需要免费注册的 `OPENALEX_API_KEY`）：用第 1 节的检索词各跑一遍（标准深度即可），结果的 BibTeX 合并到 `refs.bib`，报告网页留在本地、不必上传。
* **被引追踪**（最有效）：查谁引用了 Ye 2010（A2）、KVT 2007（B9）、Ferrando 2005（C15）、Alfassi 2007 边界力（B7）、Shin 2004（F26），重点看 2020–2026。
* 区分三类“离散对称”：样品边界（我们）、晶格 / 波导阵列（Ferrando、Kartashov 2023）、势阱各向异性（BEC）。只有第一类与我们直接竞争，后两类是要引用和对比的背景。

## 5. 本文的关键数字（核对文献时对照用）

* 对称性传递：数值 / 理论 0.99997–1.004；保角半径普适性 U(b) 差 ≤ 3.7e-4（ρ=8）。
* 圆盘 m=2 低功率窗口 b_c R² = 31.8（L=−2 分裂模）；三角形 m=2 窗口上端连续极限 2.95–2.99，Hamilton–Hopf 预言 2.93。
* 双涡旋：r_s ∝ Θ^{1/2}（渐近斜率 0.500），Θ_c = 1.176；关断后自由共转 0.1278 对预言 0.1262；间距放大偏差 ≤ 4%；
  锁定区间（第五轮 E5e，U=15）：Θ=0.125 时 −0.148…+0.029，Θ=1.1 时 −0.043…+0.007。
