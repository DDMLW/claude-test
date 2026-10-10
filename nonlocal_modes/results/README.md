# 结果上传规范（师弟必读）

## 1. 目录约定
运行任何实验脚本后会自动生成两个目录：

| 目录 | 内容 | 是否上传 |
|---|---|---|
| `nonlocal_modes/results/<实验>/<tag>/` | 图 (PNG)、表 (CSV)、`log.txt`（拆分并行时为 `log_<任务名>.txt`）、`run_info.json` | **上传** |
| `nonlocal_modes/data/<实验>/<tag>/` | 每个解的双精度原始数据 (`.npz`)，体积大 | **不上传**（已在 .gitignore 中；但请保留在本机，E3 要读 E2 的数据） |

* `<tag>`：默认等于精度档位名——正式运行 `main`，高精度 `--hq` 为 `hq`，`--quick` 冒烟测试为 `quick`。
  **上传 `main` 和 `hq`**（`quick*` 目录已被 .gitignore 忽略）。
* 超算批量运行（`experiments/hpc/`）还会生成 `hq_coarse`（粗一档网格的收敛对照）和 `hq_noisy`（E4 输入加噪声），都要上传。
* 若同一实验需要其他正式运行（例如复核某个点），请加 `--tag 名字`，各自独立成目录，都上传。
* 论文数据以 `hq` 为准；`main` 可作为补充或在没有超算时使用。

## 2. 每个实验必须上传的文件
见各实验目录下的 README.md 的“上传清单”一节。通用要求：
1. 整个 `results/<实验>/<tag>/` 目录（包括 `log*.txt` 和 `run_info*.json`，它们记录了参数、代码版本和运行时间，缺了就无法复现）；
2. 在仓库根目录 `Progress.md` 的“运行记录”表格中追加一行（日期、实验、tag、机器、耗时、是否有报错、一句话观察）。

单个结果目录一般 < 20 MB。如果某个目录超过 50 MB，请先在 Progress.md 里说明，不要直接上传大文件。

## 3. 上传方法
**方法 A：命令行（推荐）**
```bash
cd <仓库根目录>
git pull                                   # 先同步, 避免冲突
git add nonlocal_modes/results/E2_transmutation_atlas/hq Progress.md
git commit -m "E2a hq run results"
git push
```
超算上跑完的全部 hq 结果可以一次提交：`git add nonlocal_modes/results Progress.md`（`data/` 和 `quick*` 会被自动忽略）。

**方法 B：GitHub 网页**：进入 `nonlocal_modes/results/<实验>/` → Add file → Upload files → 把本机 `hq`（或 `main`）文件夹整个拖进去 → Commit。

上传到哪个分支：与代码所在分支相同（目前是 `claude/nice-planck-5809ej`；若导师已合并到 `main`，就用 `main`）。

## 4. 出问题怎么办
* 任何脚本报错：**不要改代码**，把 `log.txt` 上传，并在 Progress.md 里写清楚命令和报错的最后几行。
* 结果与 README 中“预期结果”明显不符：照常上传，并在 Progress.md 里注明“与预期不符：……”。这可能是重要发现，也可能是 bug，需要分析。
