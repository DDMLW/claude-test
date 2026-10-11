"""
run_all.py —— 按顺序运行全部实验脚本, 最后打印成功/失败汇总。

  python run_all.py --quick          先跑一遍所有 quick 模式 (约 1 小时), 确认程序在你的电脑上都能运行
  python run_all.py                  正式运行 (很久! 建议按 README 逐个实验跑, 或多开终端并行)
  python run_all.py --only E0 E1     只跑指定实验

E3 依赖 E2a 的数据 (相同 tag), 所以顺序不能乱。
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    ("E0", "E0_validation/e0a_convergence.py"),
    ("E0", "E0_validation/e0b_twin_vortex.py"),
    ("E0", "E0_validation/e0c_surface_solitons.py"),
    ("E0", "E0_validation/e0d_grid_artifacts.py"),
    ("E1", "E1_symmetry_transfer/e1a_index_anisotropy.py"),
    ("E1", "E1_symmetry_transfer/e1b_protected_modes.py"),
    ("E2", "E2_transmutation_atlas/e2a_atlas.py"),
    ("E2", "E2_transmutation_atlas/e2b_scaling.py"),
    ("E3", "E3_stability/e3_stability.py"),
    ("E4", "E4_dynamic_transmutation/e4_lg_launch.py"),
    ("E5", "E5_thermal_reconfiguration/e5a_boundary_heating.py"),
    ("E5", "E5_thermal_reconfiguration/e5b_dynamic_control.py"),
    # 第二轮 (2026-10): 正式运行只在超算上做 (hpc/make_jobs.py --phase 2); 这里只用于 --quick 冒烟测试
    ("E3", "E3_stability/e3b_instability_windows.py"),
    ("E4", "E4_dynamic_transmutation/e4b_stable_launch.py"),
    ("E5", "E5_thermal_reconfiguration/e5c_heated_stability.py"),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true")
    p.add_argument("--only", nargs="+", default=None, help="例如 E0 E2")
    args = p.parse_args()
    summary = []
    for exp, script in SCRIPTS:
        if args.only and exp not in args.only:
            continue
        cmd = [sys.executable, os.path.join(HERE, script)] + (["--quick"] if args.quick else [])
        print(f"\n######## {' '.join(cmd)}", flush=True)
        t0 = time.time()
        ret = subprocess.call(cmd, cwd=os.path.dirname(os.path.join(HERE, script)))
        summary.append((script, ret, time.time() - t0))
    print("\n================ summary ================")
    for script, ret, dt in summary:
        print(f"{'OK  ' if ret == 0 else 'FAIL'}  {dt / 60:7.1f} min  {script}")
    if any(r != 0 for _, r, _ in summary):
        print("\n有脚本失败: 请把对应的 results/<实验>/<tag>/log.txt 一并上传, 并在 Progress.md 中记录。")


if __name__ == "__main__":
    main()
