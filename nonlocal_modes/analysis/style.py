"""
style.py —— 分析/论文图的统一样式 (matplotlib)
=============================================
* 分类配色: 固定顺序 (蓝、橙、青、黄、品红、绿、紫、红), 已用色盲安全校验 (相邻色 CVD Delta E >= 9),
  不循环使用; 第 9 个系列不再生成新颜色 (折叠或分面)。
* 浅色的青、黄、品红对白底对比度 < 3:1, 所以每个系列同时用不同的标记形状 + 图例 (颜色从不是唯一的身份线索)。
* 线宽 2 px, 标记 >= 8 px, 网格为细的实线浅灰; 文字一律用中性墨色 (不用系列颜色)。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e3e2de"
SURFACE = "#fcfcfb"

# 样品的固定颜色/标记 (颜色跟随样品, 不随图中系列数目改变)
SAMPLE_ORDER = ["rect", "poly3", "poly4", "poly5", "poly6", "disk"]
SAMPLE_LABEL = {"rect": "rectangle (n=2)", "poly3": "triangle (n=3)", "poly4": "square (n=4)",
                "poly5": "pentagon (n=5)", "poly6": "hexagon (n=6)", "disk": "disk (control)"}


def sample_style(kind):
    i = SAMPLE_ORDER.index(kind)
    return dict(color=PALETTE[i], marker=MARKERS[i])


def series_style(i):
    if i >= len(PALETTE):
        raise ValueError("more than 8 series: fold into 'other' or use small multiples")
    return dict(color=PALETTE[i], marker=MARKERS[i])


def apply():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK_SECONDARY, "axes.labelcolor": INK, "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY, "text.color": INK, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "grid.linestyle": "-", "axes.axisbelow": True,
        "lines.linewidth": 1.6, "lines.markersize": 5.5, "lines.markeredgewidth": 1.0,
        "lines.markeredgecolor": SURFACE, "legend.frameon": False, "legend.fontsize": 7.5,
        "font.size": 9, "axes.titlesize": 9.5, "axes.spines.top": False, "axes.spines.right": False,
    })


def save(fig, path, dpi=160, tight=True):
    """tight=False: 版面由调用者用 gridspec/subplots_adjust 安排 (如带相位色条的场图), 保存时按内容裁边。"""
    if tight:
        fig.tight_layout()
        fig.savefig(path, dpi=dpi)
    else:
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
