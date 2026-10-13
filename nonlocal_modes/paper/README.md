# 论文正文草稿

| 文件 | 内容 |
|---|---|
| `manuscript.md` | **正文唯一源文件**，v0.2（2026-10-13，第六轮之后）：摘要、亮点、关键词、引言、第 2–9 节、附录 A（数值方法）/ B（对照与收敛）、Fig.1–6 与 SI Fig. S1–S2 图注、后置节占位 |
| `refs.bib` | 师弟第六轮核对的 BibTeX（61 条），键名即正文的 `[@key]` |
| `build.py` | 生成 Chaos, Solitons & Fractals（Elsevier elsarticle）版：pandoc → `latex/main.tex` + `body_*.tex`，pdflatex + bibtex → `latex/main.pdf` |
| `latex/` | 生成物（`main.tex`、`body_*.tex`、`refs.bib` 副本、`main.pdf`）；**不要手改**，改 `manuscript.md` 后重新 `python build.py` |

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
| Fig. S1 圆盘窗口的径向独立复核 | `results/E3_stability/hq/e3c_radial_vs_2d.png` | `experiments/E3_stability/e3c_radial_disk.py --hq` |
| Fig. S2 种子传播 | `twin_vortex_seed.png` | `a7_twin_vortex.py` |
| 其他备选 SI | `twin_vortex_SI.png`、`heated_window.png`、`triangle_convergence.png`、`stability_map.png`、`e4_growth_vs_e3.png` | `a7`、`a8`、`a5`、`a2`、`a3` |

**生成投稿版**（CSF / elsarticle；需要 pandoc，编译需要 pdflatex + bibtex + elsarticle）：
```bash
cd nonlocal_modes/paper
python build.py              # → latex/main.pdf（preprint 12pt、行号、数字引用 elsarticle-num；缺图时用占位框，不中断编译）
python build.py --no_pdf     # 只生成 tex
pandoc manuscript.md --citeproc --bibliography refs.bib -o manuscript.docx   # 给导师批注用的 Word 版
```
云端容器默认没有 TeX，每个新会话需先安装（约 1 分钟）：
`apt-get install -y --no-install-recommends texlive-latex-base texlive-latex-recommended texlive-publishers texlive-fonts-recommended`；
本地可用 TeX Live / MiKTeX，或把 `latex/` 目录与 `results/analysis/hq/`、`results/E3_stability/hq/` 中的图一起上传 Overleaf（图路径见 `main.tex` 的 `\graphicspath`）。

**待办**：
* 作者、单位、CRediT、致谢、数据可用性的仓库地址（正文后置节的占位）；
* 讨论中的材料常数（铅玻璃 dn/dT、k 是否含 n₀）请导师确认；KVT 2007 的引用措辞请导师把关（§5，Progress.md §4-13）；
* Kapitula 2013 的期刊信息有一处不一致待师弟按 DOI 核对（Progress.md §5-24）；
* 投稿信与图形摘要（CSF 可选）。
