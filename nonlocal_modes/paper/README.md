# 论文正文草稿

| 文件 | 内容 |
|---|---|
| `manuscript.md` | 英文正文草稿 v0.1（2026-10-13，第五轮分析之后）：摘要、引言、第 2–7 节结果、讨论、结论、附录 A（数值方法）/ B（对照与收敛）、Fig.1–6 图注、引用键清单 |
| `refs.bib` | （待师弟按 `docs/literature_task.md` 交付）核对过的 BibTeX，键名与正文的 `[@key]` 一致 |

**数字的来源**：正文中的每个数字都来自 `Progress.md` §3 的结论 C1–C13（括号 `[C#]` 标注出处，定稿时删除），不在正文中引入新数字。
结论有更新时先改 Progress.md，再改正文。

**图**（全部由 `nonlocal_modes/analysis/` 的脚本从已上传的 CSV 生成，在 `results/analysis/hq/`）：

| 图 | 文件 | 脚本 |
|---|---|---|
| Fig.1 对称性传递 | `symmetry_transfer.png` | `a0_symmetry_transfer.py` |
| Fig.2 嬗变图谱（场图） | `atlas_fields.png` | `a9_atlas_fields.py --tag hq --b 1.0` |
| Fig.3 分数标度律 | `scaling_local_slopes.png` | `a1_scaling_laws.py` |
| Fig.4 失稳窗口机理 | `window_mechanism.png` | `a4_window_mechanism.py` |
| Fig.5 动态嬗变与可观测窗口 | `dynamic_window.png` | `a6_dynamic_window.py` |
| Fig.6 可重构双涡旋 | `twin_vortex_fields.png`（顶行）+ `twin_vortex.png` | `a7_twin_vortex.py` |
| SI | `twin_vortex_SI.png`、`heated_window.png`、`triangle_convergence.png`、`stability_map.png`、`e4_growth_vs_e3.png` | `a7`、`a8`、`a5`、`a2`、`a3` |

**生成 LaTeX / Word**（本地或云端都有 pandoc；云端没有 TeX，不能编译 PDF）：
```bash
cd nonlocal_modes/paper
pandoc manuscript.md --standalone -o manuscript.tex          # 期刊确定后换成期刊模板 (elsarticle / revtex)
pandoc manuscript.md -o manuscript.docx                       # 给导师批注用
# refs.bib 就绪后: pandoc manuscript.md --citeproc --bibliography refs.bib -o manuscript.docx
```

**待办**（正文中以 HTML 注释标出）：
* `TODO(lit)`：引言与讨论中与先行工作的定位、KVT 2007 的核对、BEC 旋转各向异性势阱 / 旋转晶格 / azimuthon 的对照——等师弟的 `docs/literature_review.md`；
* `TODO(p6)`：第 7.3 节的种子传播检验换成第六轮 `hq_seed` 的 15 个点（云端 h=0.2 预检的数字只作占位）；
* 讨论中的材料常数（铅玻璃 dn/dT、k 是否含 n₀）请导师确认；
* 期刊与模板：按师弟的期刊建议与导师意见决定（Progress.md §7 目标期刊：中科院一区或以上）；
* 作者、单位、致谢、数据与代码可用性声明（代码在本仓库；论文数据 = `nonlocal_modes/results/`）。
