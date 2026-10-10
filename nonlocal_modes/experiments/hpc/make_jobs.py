"""
make_jobs.py —— 生成超算批量任务清单 (每行一条独立命令, 每条命令单核运行)
==========================================================================
全部实验拆成相互独立的小任务, 写成四个清单 (hpc/jobs/*.txt):

  stage1.txt  E2a 解族: 按 (样品 x 拓扑荷) 拆成 36 个任务 + 网格收敛对照 (粗一档网格 --h_scale 1.5, tag <tag>_coarse) 8 个
  stage2.txt  依赖 stage1: E2a 汇总图谱 (--plot_only); E3 稳定性 (按 样品 x 拓扑荷 拆, 读取 stage1 的解族) + 粗网格对照
  indep.txt   与上面无关的实验, 提交后立即开始: E0a-E0d, E1a, E1b, E2b, E4 (按样品拆, 另有加噪声组 <tag>_noisy), E5a, E5b
  final.txt   依赖 stage2 与 indep: E3 / E4 汇总表 (--table_only); compare_tags.py 比较 <tag> 与 <tag>_coarse (网格收敛性)

  依赖关系:  stage1 -> stage2 -> final <- indep        (最长的一条链约为 E2a 2 h + E3 8 h)

用法 (在 nonlocal_modes/experiments 目录下):
  python hpc/make_jobs.py                       # 默认 --profile hq (论文数据)
  python hpc/make_jobs.py --profile main        # 普通精度
  python hpc/make_jobs.py --python /path/to/envs/nl/bin/python   # 指定 Python 解释器
  python hpc/make_jobs.py --no_convergence      # 不做粗网格对照

然后:
  SLURM 集群:   bash hpc/submit_slurm.sh       (自动按依赖关系提交四个作业数组)
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
    lists = {"stage1": [], "stage2": [], "indep": [], "final": []}
    s1, s2, ind, fin = lists["stage1"], lists["stage2"], lists["indep"], lists["final"]
    # ---------------- stage1: E2a 解族 (E3 的输入)
    for kind in SAMPLES:
        for m in MS:
            s1.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--samples", kind, "--ms", m,
                          "--job", f"{kind}_m{m}"))
    if not a.no_convergence:
        for kind, m in CONVERGENCE_FAMILIES:
            s1.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--samples", kind, "--ms", m, *cf,
                          "--job", f"{kind}_m{m}", t=coarse))
    # ---------------- stage2: 依赖 stage1
    s2.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--plot_only", "--job", "assemble"))
    for kind in SAMPLES:
        for m in MS:
            s2.append(cmd("E3_stability/e3_stability.py", "--samples", kind, "--ms", m, "--job", f"{kind}_m{m}"))
    if not a.no_convergence:
        s2.append(cmd("E2_transmutation_atlas/e2a_atlas.py", "--plot_only", *cf, "--job", "assemble", t=coarse))
        for kind, m in CONVERGENCE_FAMILIES:
            s2.append(cmd("E3_stability/e3_stability.py", "--samples", kind, "--ms", m, *cf,
                          "--job", f"{kind}_m{m}", t=coarse))
    # ---------------- indep: 独立实验 (最长的放前面, 先开始)
    # 同一实验目录下的几个脚本同时运行, 用 --job 让它们各写各的 log_<job>.txt / run_info_<job>.json
    ind += [cmd("E0_validation/e0b_twin_vortex.py", "--job", "e0b"),
            cmd("E0_validation/e0c_surface_solitons.py", "--job", "e0c"),
            cmd("E2_transmutation_atlas/e2b_scaling.py", "--job", "e2b"),
            cmd("E5_thermal_reconfiguration/e5a_boundary_heating.py", "--job", "e5a"),
            cmd("E5_thermal_reconfiguration/e5b_dynamic_control.py", "--job", "e5b"),
            cmd("E1_symmetry_transfer/e1b_protected_modes.py", "--job", "e1b"),
            cmd("E0_validation/e0a_convergence.py", "--job", "e0a"),
            cmd("E0_validation/e0d_grid_artifacts.py", "--job", "e0d"),
            cmd("E1_symmetry_transfer/e1a_index_anisotropy.py", "--job", "e1a")]
    for kind in E4_SAMPLES:
        ind.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--samples", kind, "--job", kind))
        ind.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--samples", kind, "--noise", 0.01,
                       "--job", kind, t=f"{tag}_noisy"))
    # ---------------- final: 依赖 stage2 与 indep
    fin.append(cmd("E3_stability/e3_stability.py", "--table_only", "--job", "assemble"))
    fin.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--table_only", "--job", "assemble"))
    fin.append(cmd("E4_dynamic_transmutation/e4_lg_launch.py", "--table_only", "--job", "assemble",
                   t=f"{tag}_noisy"))
    if not a.no_convergence:
        fin.append(cmd("E3_stability/e3_stability.py", "--table_only", *cf, "--job", "assemble", t=coarse))
        fin.append(f"{py} compare_tags.py E2_transmutation_atlas {tag} {coarse}")
        fin.append(f"{py} compare_tags.py E3_stability {tag} {coarse}")

    os.makedirs(os.path.join(HERE, "jobs"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    for name in ("stage1", "stage2", "indep", "final"):
        path = os.path.join(HERE, "jobs", f"{name}.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lists[name]) + "\n")
        print(f"{name:7s}: {len(lists[name]):3d} jobs -> {os.path.relpath(path, os.path.dirname(HERE))}")
    for old in ("stage3.txt",):                          # 旧版本的清单, 避免误用
        if os.path.exists(os.path.join(HERE, "jobs", old)):
            os.remove(os.path.join(HERE, "jobs", old))
    print(f"\nprofile={a.profile}, tag={tag}" + ("" if a.no_convergence else f", grid-convergence tag={coarse}"))
    print("next: bash hpc/submit_slurm.sh   (SLURM)   or   python hpc/run_local.py -n <cores>   (no scheduler)")


if __name__ == "__main__":
    main()
