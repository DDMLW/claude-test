"""
make_jobs.py —— 生成超算批量任务清单 (每行一条独立命令, 每条命令单核运行)
==========================================================================
全部实验拆成相互独立的小任务, 写成四个清单 (hpc/jobs/*.txt):

  stage1.txt  E2a 解族: 按 (样品 x 拓扑荷) 拆成 36 个任务 + 网格收敛对照 (粗一档网格 --h_scale 1.5, tag <tag>_coarse) 8 个
  stage2.txt  依赖 stage1: E2a 汇总图谱 (--plot_only); E3 稳定性 (按 样品 x 拓扑荷 拆, 读取 stage1 的解族) + 粗网格对照
  indep.txt   与上面无关的实验, 提交后立即开始: E0a-E0d, E1a, E1b, E2b, E4 (按样品拆, 另有加噪声组 <tag>_noisy), E5a, E5b
  final.txt   依赖 stage2 与 indep: E3 / E4 汇总表 (--table_only); compare_tags.py 比较 <tag> 与 <tag>_coarse (网格收敛性)

  依赖关系:  stage1 -> stage2 -> final <- indep        (最长的一条链约为 E2a 2 h + E3 8 h)

第二轮 (--phase 2, 2026-10 hq 结果分析后的补充实验, 见 docs/research_plan.md 6.1 节) 另写两个清单:
  p2_main.txt   全部互相独立: E3b 失稳窗口 (8 个情形 + poly3:2 阈值的网格加密 h=0.15/0.075/0.05 + 旋转检验
                + 圆盘 m=2 阈值加密 + 圆盘通道谱 2 个), E5c 加热态稳定性 (3 个情形 x 8 个 Theta),
                E4b 稳定区/种子扫描/加热样品发射 (14 个 run + 2 个网格对照)
  p2_final.txt  依赖 p2_main: 各 tag 的汇总表 (--table_only)
  第二轮不重算第一轮的任何东西, 也不覆盖第一轮的清单与状态文件 (status_stage1.tsv 等)。
第三轮 (--phase 3, 第二轮结果分析后的补漏与新实验, 见 Progress.md §8) 写 p3_main.txt / p3_final.txt:
  补漏: E3b rect:2 (解族起步修正后重跑)、E5c (3:1) 全部 Theta (b 范围扩到 3.0);
  新实验: 圆盘 m=2 通道谱延伸到 b=3.4 (检验三角形宽窗口的 Krein 机理)、E5c 四极加热 (2:2) 与 C_3 加热高 b 端、
         小 Theta 线性律检验 (3:1, Theta=0.02-0.16)、E5b 四极加热双涡旋的关断/旋转、E4b 定态发射对照与四极加热发射。

用法 (在 nonlocal_modes/experiments 目录下):
  python hpc/make_jobs.py                       # 默认 --profile hq (论文数据)
  python hpc/make_jobs.py --phase 2             # 第二轮: 生成 p2_main.txt / p2_final.txt
  python hpc/make_jobs.py --phase 3             # 第三轮: 生成 p3_main.txt / p3_final.txt
  python hpc/make_jobs.py --profile main        # 普通精度
  python hpc/make_jobs.py --python /path/to/envs/nl/bin/python   # 指定 Python 解释器
  python hpc/make_jobs.py --no_convergence      # 不做粗网格对照

然后:
  SLURM 集群:   bash hpc/submit_slurm.sh       (自动按依赖关系提交四个作业数组; 第二轮: bash hpc/submit_slurm.sh p2)
  PBS 集群:     见 hpc/README.md
  无作业调度的多核服务器:  python hpc/run_local.py -n 32          (第二轮: 加 --phase 2)
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
# ---- 第二轮 (--phase 2)
E3B_CASES = ["poly3:2", "poly3:1", "disk:2", "disk:1", "poly4:2", "rect:2", "poly5:2", "poly6:2"]   # 慢的在前
# 三角形 m=2 高 b 端的弱失稳 (kappa=0, Im≈0.7) 增长率与 h=0.1 的离散误差同量级 (云端: b=2.8 处 h=0.1/0.075/0.05 为
# 稳定/2.1e-3/4.3e-3, 按 h^2 外推约 6e-3), 阈值不收敛 -> 在同一组 b 点 (2.2, 2.6, ..., 5.0) 上算四档网格, 按 h^2 外推增长率
E3B_REFINE = [  # (情形, h_scale, tag 后缀, b 区间, 额外参数): 网格加密与旋转检验, 不做非线性传播
    ("poly3:2", 1.5, "coarse", (2.2, 5.0), ["--spec_step", 0.4]),
    ("poly3:2", 0.75, "h0075", (2.2, 5.0), ["--spec_step", 0.4]),
    ("poly3:2", 0.5, "h005", (2.2, 5.0), ["--spec_step", 0.4]),   # h=0.05 每个谱点约 50 分钟, 共约 11 h
    ("poly3:2", 1.0, "", (1.8, 3.4), ["--theta0", 15]),           # 旋转检验: 多边形相对网格转 15 度 (同一 tag)
    ("disk:2", 0.5, "h005", (0.3, 0.8), []),
]
E5C_CASES = ["3:2", "4:2", "3:1"]
E5C_THETAS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
E4B_CASES = ["poly3:2:15", "poly3:2:15:1e-7", "poly3:2:15:1e-5", "poly3:2:15:1e-3",       # 与 e4b 的 ALL_CASES 相同
             "poly4:3:20:1e-14", "poly4:3:20:1e-11", "poly4:3:20:1e-8", "poly4:3:20:1e-5",
             "poly3:2:25", "poly3:2:30", "poly3:2:25:1e-2", "disk:3:20", "heat3:2:15", "heat3:2:15:1e-2"]
E4B_GRID = [("disk:3:20", 1.5, "coarse"), ("poly3:2:25", 0.5, "h005")]   # 网格对照 (圆盘破裂距离; 亚网格 r_s)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile", choices=["hq", "main"], default="hq")
    p.add_argument("--tag", default=None, help="结果标签 (默认与 profile 同名)")
    p.add_argument("--python", default="python", help="Python 解释器 (默认 python)")
    p.add_argument("--no_convergence", action="store_true", help="不生成粗网格对照任务")
    p.add_argument("--coarse_factor", type=float, default=1.5, help="对照网格 = 本档网格 x 此因子 (默认 1.5)")
    p.add_argument("--phase", type=int, choices=[1, 2, 3], default=1,
                   help="1: 第一轮全部实验 (stage1/stage2/indep/final); 2: 第二轮 (p2_main/p2_final); "
                        "3: 第三轮 (p3_main/p3_final)")
    a = p.parse_args()
    tag = a.tag or a.profile
    coarse = f"{tag}_coarse"
    flag = "--hq" if a.profile == "hq" else ""
    py = a.python

    def cmd(script, *extra, t=tag):
        parts = [py, script] + ([flag] if flag else []) + ["--tag", t] + [str(x) for x in extra]
        return " ".join(parts)

    if a.phase in (2, 3):
        lists = phase2_lists(cmd, tag, py) if a.phase == 2 else phase3_lists(cmd, tag)
        write_lists(lists, tuple(lists))
        print(f"\nphase {a.phase}, profile={a.profile}, tag={tag}")
        print(f"next: bash hpc/submit_slurm.sh p{a.phase}   (SLURM)   or   "
              f"python hpc/run_local.py -n <cores> --phase {a.phase}")
        return
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

    write_lists(lists, ("stage1", "stage2", "indep", "final"))
    for old in ("stage3.txt",):                          # 旧版本的清单, 避免误用
        if os.path.exists(os.path.join(HERE, "jobs", old)):
            os.remove(os.path.join(HERE, "jobs", old))
    print(f"\nprofile={a.profile}, tag={tag}" + ("" if a.no_convergence else f", grid-convergence tag={coarse}"))
    print("next: bash hpc/submit_slurm.sh   (SLURM)   or   python hpc/run_local.py -n <cores>   (no scheduler)")


def write_lists(lists, names):
    os.makedirs(os.path.join(HERE, "jobs"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    for name in names:
        path = os.path.join(HERE, "jobs", f"{name}.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lists[name]) + "\n")
        print(f"{name:8s}: {len(lists[name]):3d} jobs -> {os.path.relpath(path, os.path.dirname(HERE))}")


def phase2_lists(cmd, tag, py):
    """第二轮: p2_main (全部独立, 慢的在前) 与 p2_final (汇总)。"""
    main, fin = [], []
    e3b, e5c, e4b = ("E3_stability/e3b_instability_windows.py", "E5_thermal_reconfiguration/e5c_heated_stability.py",
                     "E4_dynamic_transmutation/e4b_stable_launch.py")
    e3b_tags, e4b_tags = {tag}, {tag}
    # E3b 网格加密 (h=0.05 最慢, 放最前)
    for case, hs, suffix, (b1, b2), extra in sorted(E3B_REFINE, key=lambda r: r[1]):
        t = f"{tag}_{suffix}" if suffix else tag
        e3b_tags.add(t)
        job = case.replace(":", "_m") + (f"_rot{extra[1]}" if "--theta0" in extra else "") + "_refine"
        main.append(cmd(e3b, "--cases", case, "--b_range", b1, b2, "--no_evolution",
                        *(["--h_scale", hs] if hs != 1.0 else []), *extra, "--job", job, t=t))
    for case in E3B_CASES:
        main.append(cmd(e3b, "--cases", case, "--job", case.replace(":", "_m")))
    for case in ("disk:2", "disk:1"):
        main.append(cmd(e3b, "--cases", case, "--channels", "--job", case.replace(":", "_m") + "_channels"))
    # E5c: 每个 (情形, Theta) 一个任务
    for case in E5C_CASES:
        n, m = case.split(":")
        for T in E5C_THETAS:
            main.append(cmd(e5c, "--cases", case, "--thetas", T, "--job", f"e5c_n{n}_m{m}_T{T:g}"))
    # E4b: 每个 run 一个任务, 另加网格对照
    for case in E4B_CASES:
        main.append(cmd(e4b, "--cases", case, "--job", "e4b_" + case.replace(":", "_")))
    for case, hs, suffix in E4B_GRID:
        e4b_tags.add(f"{tag}_{suffix}")
        main.append(cmd(e4b, "--cases", case, "--h_scale", hs, "--job", "e4b_" + case.replace(":", "_"),
                        t=f"{tag}_{suffix}"))
    # 汇总
    for t in sorted(e3b_tags):
        fin.append(cmd(e3b, "--table_only", "--job", "assemble", t=t))
    fin.append(cmd(e5c, "--table_only", "--job", "assemble"))
    for t in sorted(e4b_tags):
        fin.append(cmd(e4b, "--table_only", "--job", "assemble", t=t))
    return {"p2_main": main, "p2_final": fin}


def phase3_lists(cmd, tag):
    """第三轮: p3_main (全部独立, 慢的在前) 与 p3_final (汇总)。"""
    main, fin = [], []
    e3b, e5c, e5b, e4b = ("E3_stability/e3b_instability_windows.py", "E5_thermal_reconfiguration/e5c_heated_stability.py",
                          "E5_thermal_reconfiguration/e5b_dynamic_control.py",
                          "E4_dynamic_transmutation/e4b_stable_launch.py")
    hib, small = f"{tag}_hib", f"{tag}_smallT"
    # E5c C_3 加热 m=2 的高 b 端: 第二轮 Theta>=1 时失稳区一直延伸到扫描上端 b=3.0, 找上端 b_up(Theta) 与稳定的加热态
    for T in (1.0, 1.5, 2.0, 2.5, 3.0):
        main.append(cmd(e5c, "--cases", "3:2", "--thetas", T, "--b_range", 2.9, 6.0, "--job", f"e5c_n3_m2_T{T:g}_hib",
                        t=hib))
    # E5c 四极加热 (双涡旋, 云端 h=0.2 检查: Theta=0.25、1 在 U=15 附近稳定) 与 (3:1) 重跑 (b 范围扩到 3.0)
    for case, thetas in (("2:2", (0.0, 0.03, 0.06, 0.125, 0.25, 0.5, 0.75, 1.0)),
                         ("3:1", (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))):
        n, m = case.split(":")
        for T in thetas:
            main.append(cmd(e5c, "--cases", case, "--thetas", T, "--job", f"e5c_n{n}_m{m}_T{T:g}"))
    # E3b: rect:2 重跑 (解族起步修正); 圆盘 m=2 通道谱延伸到 b=3.4 (另一个 tag, 不覆盖 b<=1.35 的结果)
    main.append(cmd(e3b, "--cases", "rect:2", "--job", "rect_m2"))
    main.append(cmd(e3b, "--cases", "disk:2", "--channels", "--b_range", 1.4, 3.4, "--job", "disk_m2_channels",
                    t=hib))
    # 小 Theta 线性律: C_3 加热 m=1, 圆盘通道 L=-2 与 +1 的交叉点 b≈0.26 (R=10) 附近, 细 b 网格 + 只扫失稳模的频段
    for T in (0.0, 0.02, 0.04, 0.08, 0.16):
        main.append(cmd(e5c, "--cases", "3:1", "--thetas", T, "--b_range", 0.1, 0.5, "--spec_step", 0.01,
                        "--omega", 0.02, 0.2, "--job", f"e5c_n3_m1_T{T:g}_small", t=small))
    # E5b: 四极加热双涡旋的关断/旋转 (U=15)
    for T in (0.5, 1.0):
        main.append(cmd(e5b, "--n", 2, "--Theta", T, "--job", f"e5b_n2_T{T:g}"))
    # E4b: 定态发射对照 (第二轮 A 组慢增长的来源) 与四极加热中的 LG 发射 (动态形成双涡旋)
    for case in ("poly3:2:25:0:stat", "poly3:2:25:1e-6:stat"):
        main.append(cmd(e4b, "--cases", case, "--job", "e4b_" + case.replace(":", "_")))
    for case in ("heat2:2:15", "heat2:2:15:1e-2"):
        main.append(cmd(e4b, "--cases", case, "--Theta", 0.5, "--job", "e4b_" + case.replace(":", "_")))
    # 汇总
    fin.append(cmd(e3b, "--table_only", "--job", "assemble"))
    for t in (tag, hib, small):
        fin.append(cmd(e5c, "--table_only", "--job", "assemble", t=t))
    fin.append(cmd(e4b, "--table_only", "--job", "assemble"))
    return {"p3_main": main, "p3_final": fin}


if __name__ == "__main__":
    main()
