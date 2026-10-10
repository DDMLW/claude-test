"""
compare_tags.py —— 比较同一实验两个 tag 的结果表 (网格收敛性检验)
================================================================
典型用法: hq 档 (h=0.1, 大样品 0.15) 与粗一档网格 (--h_scale 1.5, tag hq_coarse) 的对比:
  python compare_tags.py E2_transmutation_atlas hq hq_coarse
  python compare_tags.py E3_stability hq hq_coarse

做法: 对两个目录中同名的每个 CSV, 按“键列” (sample, n, m, b, rho, Theta, ... 中存在的那些) 配对各行,
逐列比较:
  * 连续量 (U, peak, n_max, ring_radius, max_growth, r_s_fit ...): 最大绝对差、最大相对差;
    二阶格式下细网格结果的误差估计 ≈ 差值 / (ratio^2 - 1) (ratio = 1.5 时为 差值 x 0.8)
  * 整数量 (central_charge, n_satellites, total_charge ...): 不一致的行数。对解族摘要 (有 r_s_fit 列) 区分
    “已分辨”的行 (不预期嬗变, 或两套网格上嬗变核半径 r_s_fit 都 >= 2h) 与“未分辨”的行 (卫星环小于网格, 网格上读出的
    中心荷/卫星数取决于 h, 不一致是正常的)。**已分辨行的拓扑量必须完全一致**。

输出: 屏幕表格 + results/<实验>/<tagA>/compare_<tagA>_vs_<tagB>.csv (随 tagA 的结果一起上传)。

判读 (写进 Progress.md):
  U、ring_radius 的相对差 < 1e-2, max_growth 的绝对差 < 1e-3 且稳定/不稳定判定一致 -> 网格已收敛;
  否则在论文中只使用 hq 结果并注明误差估计, 或在该处再加密一档 (--h 0.075) 复核。
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from nlmodes import io_utils  # noqa: E402

KEY_COLUMNS = ["sample", "case", "test", "name", "n", "m", "b", "rho", "Theta", "dz", "stencil"]
# 由网格检测的拓扑量: 嬗变核小于网格时其读数依赖 h, 只在“已分辨”的行上要求一致
GRID_TOPOLOGY = {"central_charge", "n_satellites", "satellite_charge", "W_small", "W_core", "unresolved_charge",
                 "central_num", "satellite_charge_num", "n_satellites_detected"}
IGNORE = {"iterations", "time", "n_eigs", "converged", "h", "residual", "n_solutions", "xc", "yc"}  # xc,yc: 对称性下为 0 (噪声)


def _num(v):
    if v in ("True", "False"):                 # 布尔列 (例如 E3 的 unstable) 按 1/0 比较
        return 1.0 if v == "True" else 0.0
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _key(row, keys):
    out = []
    for k in keys:
        x = _num(row.get(k))
        out.append(f"{x:.4f}" if x is not None else str(row.get(k)))
    return tuple(out)


def _fmt_key(row, keys):
    out = []
    for k in keys:
        x = _num(row.get(k))
        out.append(f"{k}={x:.4g}" if x is not None and k not in ("sample", "case", "test", "name")
                   else f"{k}={row.get(k)}")
    return ";".join(out)


def _resolved(x, y):
    """该行 (两套网格) 的嬗变结构是否已被网格分辨; 缺少所需列时返回 None (不区分)。"""
    try:
        m, v = int(float(x["m"])), int(float(x["predicted_central"]))
        if m == v:
            return True                      # 不预期嬗变 (中心荷 = m)
        return all(float(r["r_s_fit"]) >= 2 * float(r["h"]) for r in (x, y))
    except (KeyError, TypeError, ValueError):
        return None


def compare_file(fa, fb, ratio):
    ra, rb = io_utils.read_csv(fa), io_utils.read_csv(fb)
    if not ra or not rb:
        return []
    cols_a, cols_b = list(ra[0].keys()), list(rb[0].keys())
    keys = [k for k in KEY_COLUMNS if k in cols_a and k in cols_b]
    if "b" not in keys and "U" in cols_a and "U" in cols_b and "zmax" in cols_a:
        keys.append("U")                      # E4: 输入能流 U 是参数而不是结果
    if not keys:
        return [dict(column="(no key columns)", n_matched=0)]
    ib = {_key(r, keys): r for r in rb}
    pairs = [(r, ib[_key(r, keys)]) for r in ra if _key(r, keys) in ib]
    out = []
    for c in cols_a:
        if c in keys or c in IGNORE or c not in cols_b:
            continue
        vals = [(_num(x[c]), _num(y[c]), x) for x, y in pairs]
        vals = [(a, b, x) for a, b, x in vals if a is not None and b is not None]
        if not vals:
            continue
        integer = all(float(a).is_integer() and float(b).is_integer() for a, b, _ in vals)
        if integer:
            pairs_ok = [(x, y) for x, y in pairs if _num(x[c]) is not None and _num(y[c]) is not None]
            bad = [(x, y) for x, y in pairs_ok if _num(x[c]) != _num(y[c])]
            res = dict(column=c, kind="integer", n_matched=len(vals), n_mismatch=len(bad))
            if c in GRID_TOPOLOGY and bad and _resolved(*bad[0]) is not None:
                bad_res = [(x, y) for x, y in bad if _resolved(x, y)]
                res.update(n_mismatch_resolved=len(bad_res), n_mismatch_unresolved=len(bad) - len(bad_res))
                bad = bad_res
            else:
                res["n_mismatch_resolved"] = len(bad)
            res["first_mismatch"] = _fmt_key(bad[0][0], keys) if bad else ""
            out.append(res)
        else:
            d = [abs(a - b) for a, b, _ in vals]
            rel = [abs(a - b) / max(abs(a), abs(b), 1e-12) for a, b, _ in vals]
            j = max(range(len(d)), key=d.__getitem__)
            out.append(dict(column=c, kind="float", n_matched=len(vals), max_abs_diff=max(d),
                            max_rel_diff=max(rel), est_error_fine=max(d) / (ratio ** 2 - 1),
                            at=_fmt_key(vals[j][2], keys)))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("experiment", help="实验目录名, 例如 E2_transmutation_atlas")
    p.add_argument("tag_a", help="细网格 tag (例如 hq)")
    p.add_argument("tag_b", help="粗网格 tag (例如 hq_coarse)")
    p.add_argument("--ratio", type=float, default=1.5, help="两套网格间距之比 h_b/h_a (默认 1.5)")
    a = p.parse_args()
    da = os.path.join(io_utils.PROJECT_ROOT, "results", a.experiment, a.tag_a)
    db = os.path.join(io_utils.PROJECT_ROOT, "results", a.experiment, a.tag_b)
    for d in (da, db):
        if not os.path.isdir(d):
            sys.exit(f"directory not found: {d}")
    files = sorted(f for f in os.listdir(da) if f.endswith(".csv") and not f.startswith("compare_")
                   and os.path.exists(os.path.join(db, f)))
    rows = []
    print(f"comparing {a.experiment}: {a.tag_a} vs {a.tag_b} ({len(files)} common CSV files)")
    for f in files:
        res = compare_file(os.path.join(da, f), os.path.join(db, f), a.ratio)
        if not res:
            continue
        print(f"\n{f}")
        for r in res:
            r = dict(file=f, **r)
            rows.append(r)
            if r.get("kind") == "integer":
                flag = "" if r["n_mismatch_resolved"] == 0 else "   <-- MISMATCH (resolved rows)"
                extra = (f" (resolved rows {r['n_mismatch_resolved']}, unresolved rows {r['n_mismatch_unresolved']}"
                         f" - expected)" if "n_mismatch_unresolved" in r else "")
                print(f"   {r['column']:24s} integer  rows {r['n_matched']:4d}  mismatches {r['n_mismatch']}{extra}"
                      f"{('  first at ' + r['first_mismatch']) if r['first_mismatch'] else ''}{flag}")
            elif r.get("kind") == "float":
                print(f"   {r['column']:24s} float    rows {r['n_matched']:4d}  max|diff| {r['max_abs_diff']:.2e}  "
                      f"max rel {r['max_rel_diff']:.2e}  (est. error of {a.tag_a}: {r['est_error_fine']:.1e})")
            else:
                print(f"   {r['column']}")
    out = os.path.join(da, f"compare_{a.tag_a}_vs_{a.tag_b}.csv")
    io_utils.write_csv(out, rows)
    print(f"\nwritten: {out}")


if __name__ == "__main__":
    main()
