# CLAUDE.md —— 给以后在本仓库工作的 Claude 的说明

本仓库当前承载一个科研项目：**非局域非线性材料中特殊模态的设计**（边界对称性主导的涡旋孤子拓扑荷嬗变）。
用户是研究生，导师是热非局域孤子方向的论文作者；用户的师弟负责在本地运行代码并上传结果。用户用中文交流。

## 每次开始工作
1. **先读根目录 `Progress.md`**——它是唯一事实来源（已确认结论、开放问题、已知陷阱、运行记录、决策日志）。不要凭记忆。
2. 需要细节时再读 `nonlocal_modes/docs/`（research_plan / theory / numerics / references）和各实验目录的 README.md。
3. `git pull` 后查看 `nonlocal_modes/results/<实验>/main/` 中新上传的 CSV/PNG/log.txt。

## 约定
* 代码在 `nonlocal_modes/`；核心包 `nlmodes/` 改动后必须 `cd nonlocal_modes && python -m pytest -q tests` 全部通过。
* 新实验脚本放 `nonlocal_modes/experiments/E*/`，用 `common.standard_args / finalize_args`（自动写 log.txt、run_info.json），
  必须有 `--quick` 模式并在提交前实际跑通；每个实验目录要有 README（运行命令、耗时、输出、预期结果、上传清单）。
* 注释与文档用中文；`print` 只用 ASCII（Windows 控制台编码问题）；图中文字用英文。
* 原始数据 `nonlocal_modes/data/` 与 quick 结果不进 git（见 .gitignore）。
* **数值结论的纪律**（见 numerics.md 第 6 节）：卫星环 < 2h 不下结构结论；n=5,6 以 E5a 为准；结论需通过 E0d 旋转检验与网格加密；
  中心多重荷只用圆周绕数判断。理论与数值不符时先查数值，再改理论，并记入 Progress.md。
* 每次有进展（结果分析、代码修改、新决定）都更新 `Progress.md`（结论表、开放问题、运行记录、决策日志、最后更新日期）。
* 开发分支：`claude/nice-planck-5809ej`；未经用户要求不要建 PR。
