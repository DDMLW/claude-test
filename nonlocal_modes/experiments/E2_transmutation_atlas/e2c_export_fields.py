"""
E2c —— 导出论文场图 (Fig.2 嬗变图谱) 所需的小文件: 由 E2a 保存在 data/ 中的解裁剪、降为单精度后写入 results/ (上传)
=====================================================================================================
E2a 把每个解族的全部解以双精度存在 data/E2_transmutation_atlas/<tag>/<样品>_m<m>_b<b>.npz (不上传, 每个约 1-5 MB)。
论文作图脚本 (analysis/a9_atlas_fields.py) 只读 results/, 所以这里对每个 (样品, m) 取 b 最接近 --b 的解,
裁剪到光束所在的中心区域 (半宽 = --zoom 倍亮环半径 + 1, 不超过样品), 存为:
    results/E2_transmutation_atlas/<tag>/fields/<样品>_m<m>_b<b>.npz
        x, y      裁剪后的格心坐标 (float32)          psi   复光场 (complex64)       n   折射率 (float32)
        b, U, h   传播常数、能流、网格间距             sample, m, nsym, ring_radius
    results/E2_transmutation_atlas/<tag>/fields/index.csv   每个文件一行 (b、U、裁剪半宽、峰值、文件大小)
不做任何求解, 几秒钟完成 (hq 的 36 个解约 3-5 MB)。

运行:
  python e2c_export_fields.py --hq                     (默认 b=1.0, 全部样品, m=1..6; 读取 data/E2_transmutation_atlas/hq)
  python e2c_export_fields.py --hq --b 0.5 1.0 2.0     (每个 (样品, m) 导出多个 b)
  python e2c_export_fields.py --quick                  (读取 quick 档的数据, 用于检查)
"""

import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import finalize_args, standard_args  # noqa: E402

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import io_utils  # noqa: E402
from nlmodes.geometry import Grid  # noqa: E402

SAMPLES = ["rect", "poly3", "poly4", "poly5", "poly6", "disk"]


def saved_grid(d):
    """由保存的 x, y, h 重建 Grid (格心均匀网格)。"""
    g = Grid.__new__(Grid)
    g.x, g.y, g.h = np.asarray(d["x"], float), np.asarray(d["y"], float), float(d["h"])
    g.X, g.Y = np.meshgrid(g.x, g.y)
    g.shape = g.X.shape
    g.R, g.TH = np.hypot(g.X, g.Y), np.arctan2(g.Y, g.X)
    return g


def main():
    p = standard_args(__doc__)
    p.add_argument("--samples", nargs="+", default=SAMPLES)
    p.add_argument("--ms", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    p.add_argument("--b", nargs="+", type=float, default=[1.0], help="导出 b 最接近这些值的解")
    p.add_argument("--zoom", type=float, default=1.8, help="裁剪半宽 = zoom x 亮环半径 + 1")
    p.add_argument("--data_tag", default=None, help="读取哪个 tag 的 data/ (默认与 --tag 相同)")
    args = finalize_args(p, "E2_transmutation_atlas")
    src = io_utils.data_dir("E2_transmutation_atlas", args.data_tag or args.tag)
    out = io_utils.ensure_dir(os.path.join(args.res_dir, "fields"))
    index = []
    for kind in args.samples:
        for m in args.ms:
            files = glob.glob(os.path.join(src, f"{kind}_m{m}_b*.npz"))
            bs = {float(re.search(r"_b([0-9.]+)\.npz$", f).group(1)): f for f in files}
            if not bs:
                print(f"  {kind} m={m}: no saved solutions in {src}")
                continue
            for bt in args.b:
                bsel = min(bs, key=lambda v: abs(v - bt))
                d = io_utils.load_solution(bs[bsel])
                g = saved_grid(d)
                q = d["q"]
                ring = dg.ring_radius(q, g, (0.0, 0.0))
                half = min(args.zoom * ring + 1.0, g.x.max(), g.y.max())
                ix = np.where(np.abs(g.x) <= half)[0]
                iy = np.where(np.abs(g.y) <= half)[0]
                sl = (slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1))
                name = f"{kind}_m{m}_b{bsel:.2f}"
                fn = os.path.join(out, name + ".npz")
                np.savez_compressed(fn, x=g.x[sl[1]].astype(np.float32), y=g.y[sl[0]].astype(np.float32),
                                    psi=q[sl].astype(np.complex64), n=np.asarray(d["n"])[sl].astype(np.float32),
                                    b=bsel, U=d["U"], h=g.h, sample=kind, m=m,
                                    nsym=int(d["meta"].get("nsym") or 0), ring_radius=ring)
                index.append(dict(file=name + ".npz", sample=kind, m=m, b=bsel, U=d["U"], h=g.h, half_width=half,
                                  ring_radius=ring, peak=float(np.abs(q).max()), kB=round(os.path.getsize(fn) / 1024, 1)))
                print(f"  {name}: crop |x|,|y| <= {half:.2f} ({len(ix)}x{len(iy)} points), {index[-1]['kB']} kB",
                      flush=True)
    io_utils.write_csv(os.path.join(out, "index.csv"), index)
    print(f"\nE2c done: {len(index)} fields, {sum(r['kB'] for r in index) / 1024:.1f} MB -> {out}")


if __name__ == "__main__":
    main()
