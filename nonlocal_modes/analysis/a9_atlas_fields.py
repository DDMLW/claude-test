"""
a9_atlas_fields.py —— 拓扑荷嬗变图谱的场图 (论文 Fig.2)
=====================================================
读取 E2c 导出的小文件 results/E2_transmutation_atlas/<tag>/fields/<样品>_m<m>_b<b>.npz (由 E2a 的解裁剪而来)
与 E2a 的 atlas_table.csv (卫星环是否可分辨), 画 样品 x 拓扑荷 的图谱:
  每格: 色相 = 相位 (HSV), 亮度 = sqrt(|psi| / max|psi|) (暗区也能看到相位结构);
        白圈 = +1、黑圈 = -1 单位奇点 (方格法, 中心 1.25h 以外); 左上角 = 中心荷 (0.5h 圆周绕数);
        卫星环 < 2h (未分辨, atlas_table 的 resolved=False) 的格子右下角标 "unresolved"。
输出: results/analysis/<tag>/atlas_fields.png
运行:  python a9_atlas_fields.py [--tag hq] [--b 1.0]      (先在超算上运行 E2c 导出场)
"""

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

import fieldplot  # noqa: E402
import style  # noqa: E402

ROOT = os.path.join(io_utils.PROJECT_ROOT, "results")
SAMPLES = ["rect", "poly3", "poly4", "poly5", "poly6", "disk"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="hq")
    p.add_argument("--b", type=float, default=1.0, help="每格取 b 最接近这个值的导出场")
    a = p.parse_args()
    src = os.path.join(ROOT, "E2_transmutation_atlas", a.tag, "fields")
    files = glob.glob(os.path.join(src, "*_m*_b*.npz"))
    if not files:
        sys.exit(f"no exported fields in {src}: run experiments/E2_transmutation_atlas/e2c_export_fields.py first")
    have = {}
    for f in files:
        mt = re.search(r"([a-z0-9]+)_m(\d+)_b([0-9.]+)\.npz$", os.path.basename(f))
        kind, m, b = mt.group(1), int(mt.group(2)), float(mt.group(3))
        if (kind, m) not in have or abs(b - a.b) < abs(have[(kind, m)][0] - a.b):
            have[(kind, m)] = (b, f)
    resolved = {}
    ft = os.path.join(ROOT, "E2_transmutation_atlas", a.tag, "atlas_table.csv")
    if os.path.exists(ft):
        for r in io_utils.read_csv(ft):
            resolved[(r["sample"], int(r["m"]))] = r.get("resolved", "True") == "True"
    kinds = [k for k in SAMPLES if any(k == kk for kk, _ in have)]
    ms = sorted({m for _, m in have})
    style.apply()
    nr, nc = len(kinds), len(ms)
    fig = style.plt.figure(figsize=(1.9 * nc + 0.6, 1.9 * nr + 0.7))
    gs = fig.add_gridspec(nr + 1, nc, height_ratios=[1] * nr + [0.06], hspace=0.12, wspace=0.06, top=1 - 0.55 / (1.9 * nr),
                          bottom=0.03)
    for i, kind in enumerate(kinds):
        for j, m in enumerate(ms):
            ax = fig.add_subplot(gs[i, j])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if i == 0:
                ax.set_title(f"m = {m}", fontsize=10)
            if j == 0:
                ax.set_ylabel(style.SAMPLE_LABEL[kind], fontsize=8.5)
            if (kind, m) not in have:
                ax.text(0.5, 0.5, "no solution\nat this b", ha="center", va="center", fontsize=7,
                        color=style.INK_SECONDARY, transform=ax.transAxes)
                continue
            g, q, meta = fieldplot.load_field(have[(kind, m)][1])
            _, central = fieldplot.draw_field(ax, g, q)
            if j == 0:
                ax.set_ylabel(style.SAMPLE_LABEL[kind], fontsize=8.5)
            ax.text(0.04, 0.96, f"{central:+d}" if central else "0", transform=ax.transAxes, ha="left", va="top",
                    fontsize=8, color="white")
            if resolved.get((kind, m)) is False:
                ax.text(0.96, 0.04, "unresolved", transform=ax.transAxes, ha="right", va="bottom", fontsize=6,
                        color="white")
    fieldplot.phase_bar(fig.add_subplot(gs[nr, (nc - 1) // 2]))
    fig.suptitle(f"hue = phase, brightness = sqrt(amplitude) at b ~ {a.b:g}\nwhite/black circles: +1/-1 "
                 f"phase singularities; corner: central charge", fontsize=8.5, y=0.995)
    out = io_utils.ensure_dir(os.path.join(ROOT, "analysis", a.tag))
    style.save(fig, os.path.join(out, "atlas_fields.png"), dpi=170, tight=False)
    print(f"written to {out}")


if __name__ == "__main__":
    main()
