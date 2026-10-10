"""
make_jobs.py —— 生成超算批量任务清单 (每行一条独立命令, 每条命令单核运行)
==========================================================================
所有实验被拆成相互独立的小任务, 分三个阶段 (后一阶段要等前一阶段全部结束):

  stage1  E0a-E0d, E1a, E1b, E2a (按 样品 x 拓扑荷 拆成 36 个任务), E2b, E4 (按样品拆), E5a, E5b
          + 网格收敛性检验: 关键解族在粗一档网格 (--h_scale 1.5, tag <tag>_coarse) 上重算
  stage2  E2a 汇总图谱 (--plot_only); E3 (按 样品 x 拓扑荷 拆, 读取 stage1 中 E2a 的解); E4 汇总表
  stage3  E3 汇总表 (--table_only); compare_tags.py 比较 <tag> 与 <tag>_coarse (网格收敛性)

用法 (在 nonlocal_modes/experiments 目录下):
  python hpc/make_jobs.py                       # 默认 --profile hq (论文数据), 生成 hpc/jobs/stage{1,2,3}.txt
  python hpc/make_jobs.py --profile main        # 普通精度
  python hpc/make_jobs.py --python /path/to/envs/nl/bin/python   # 指定 Python 解释器
  python hpc/make_jobs.py --no_convergence      # 不做粗网格对照

然后:
  SLURM 集群:   bash hpc/submit_slurm.sh       (自动按阶段依赖提交三个作业数组)
  PBS 集群:     见 hpc/README.md
  无作业调度的多核服务器:  python hpc/run_local.py -n 32
"""

import argparse
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = ["rect", "poly3", "poly4", "poly5", "poly6", "disk"]
MS = [1, 2, 3, 4, 5, 6]
E4_SAMPLES = ["poly3", "poly4", "poly5", "poly6", "disk"]
# 网格收敛性检验所用的关键解族 (样品, m): 嬗变的典型情形 + 双涡旋 + 圆盘对照
CONVERGENCE_FAMILIES = [("poly3", 2), ("poly3", 3), ("poly4", 2), ("poly4", 3), ("poly5", 2),
                        ("rect", 2), ("rect", 3), ("disk", 2)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile", choices=["hq", "main"], default="hq")
    p.add_argument("--tag", default=None, help="结果标签 (默认与 profile 同名)")
    p.add_argument("--python", default="python", help="Python 解释器 (默认 python)")
    p.add_argument("--no_convergence", action="store_true", help="不生成粗网格对照任务")
    p.add_argument("--coarse_factor", type=float, default=1.5, help="对照网格 = 本档网格 x 此因子 (默认 1.5)")
    a = p.parse_args()
    tag = a.tag or a.profile
    coarse = f"{tag}_coarse"
    flag = "--hq" if a.profile == "hq" else ""
    py = a.python

    def cmd(script, *extra, t=tag):
        parts = [py, script] + ([flag] if flag else []) + ["--tag", t] + [str(x) for x in extra]
        return " ".join(parts)

    cf = ["--h_scale", a.coarse_factor]
    s1, s2, s3 = [], [], []
    # ---------------- stage 1
    # 同一实验目录下的几个脚本同时运行, 用 --job 让它们各写各的 log_<job>.txt / run_info_<job>.json
    s1 += [cmd("E0_validation/e0a_convergence.py", "--job", "e0a"),
           cmd("E0_validation/e0b_twin_vortex.py", "--job", "e0b"),
           cmd("E0_validation/e0c_surface_solitons.py", "--job", "e0c"),
           cmd("E0_validation/e0d_grid_artifacts.py", "--job", "e0d"),
           cmd("E1_symmetry_transfer/e1a_index_anisotropy.py", "--job", "e1a"),
           cmd("E1_symmetry_transfer/e1b_protected_modes.py", "--job", "e1b"),
           cmd("E2_transmutation_atlas/e2b_scaling.py", "--job", "e2b"),
           cmd("E5_thermal_reconfiguration/e5a_boundary_heating.py", "--job", "e5a"),
           cmd("E5_thermal_reconfiguration/e5b_dynamic_control.py", "--job", "e5b")]
    for kind in SAMPLES:
        for m in MS:
            s1.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--samples", kind, "--ms", m,
                          "--job", f"{kind}_m{m}"))
    for kind in E4_SAMPLES:
        s1.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--samples", kind, "--job", kind))
        s1.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--samples", kind, "--noise", 0.01,
                      "--job", kind, t=f"{tag}_noisy"))
    if not a.no_convergence:
        for kind, m in CONVERGENCE_FAMILIES:
            s1.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--samples", kind, "--ms", m, *cf,
                          "--job", f"{kind}_m{m}", t=coarse))
    # ---------------- stage 2
    s2.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--plot_only", "--job", "assemble"))
    for kind in SAMPLES:
        for m in MS:
            s2.append(cmd("E3_stability/e3_stability.py", "--samples", kind, "--ms", m, "--job", f"{kind}_m{m}"))
    s2.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--table_only", "--job", "assemble"))
    s2.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--table_only", "--job", "assemble",
                  t=f"{tag}_noisy"))
    if not a.no_convergence:
        s2.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--plot_only", *cf, "--job", "assemble", t=coarse))
        for kind, m in CONVERGENCE_FAMILIES:
            s2.append(cmd("E3_stability/e3_stability.py", "--samples", kind, "--ms", m, *cf,
                          "--job", f"{kind}_m{m}", t=coarse))
    # ---------------- stage 3
    s3.append(cmd("E3_stability/e3_stability.py", "--table_only", "--job", "assemble"))
    if not a.no_convergence:
        s3.append(cmd("E3_stability/e3_stability.py", "--table_only", *cf, "--job", "assemble", t=coarse))
        s3.append(f"{py} compare_tags.py E2_transmutation_atlas {tag} {coarse}")
        s3.append(f"{py} compare_tags.py E3_stability {tag} {coarse}")

    os.makedirs(os.path.join(HERE, "jobs"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    for i, lines in enumerate((s1, s2, s3), start=1):
        path = os.path.join(HERE, "jobs", f"stage{i}.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        print(f"stage{i}: {len(lines):3d} jobs -> {os.path.relpath(path, os.path.dirname(HERE))}")
    print(f"\nprofile={a.profile}, tag={tag}" + ("" if a.no_convergence else f", grid-convergence tag={coarse}"))
    print("next: bash hpc/submit_slurm.sh   (SLURM)   or   python hpc/run_local.py -n <cores>   (no scheduler)")


if __name__ == "__main__":
    main()
