"""
单元测试 (pytest)。运行:  cd nonlocal_modes && python -m pytest -q tests
全部通过才说明代码在你的电脑上工作正常; 任何一个失败都要先停下来报告, 不要继续跑实验。
总耗时约 1-2 分钟。
"""

import os
import sys

import numpy as np
import pytest
from scipy.special import exp1

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nlmodes import diagnostics as dg  # noqa: E402
from nlmodes import guesses, newton, propagate, stability, theory  # noqa: E402
from nlmodes.geometry import (disk_sample, polygon_sample, rectangle_sample,  # noqa: E402
                              surface_interface_sample, harmonic_multipole)
from nlmodes.model import NonlocalModel  # noqa: E402
from nlmodes.operators import laplacian  # noqa: E402


# ---------------------------------------------------------------------------
# 1. 差分算子: 与离散本征值的精确公式比较 (机器精度)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("stencil", [5, 9])
def test_laplacian_spectrum_rectangle(stencil):
    h, nx, ny = 0.5, 7, 5
    mask = np.zeros((ny + 2, nx + 2), bool)
    mask[1:-1, 1:-1] = True
    D = laplacian(mask, h, stencil).toarray()
    ev = np.sort(np.linalg.eigvalsh(D))
    a = np.pi * np.arange(1, nx + 1) / (nx + 1)
    b = np.pi * np.arange(1, ny + 1) / (ny + 1)
    A, B = np.meshgrid(a, b)
    if stencil == 5:
        ex = -(4 / h ** 2) * (np.sin(A / 2) ** 2 + np.sin(B / 2) ** 2)
    else:
        ex = (4 * (2 * np.cos(A) + 2 * np.cos(B)) + 4 * np.cos(A) * np.cos(B) - 20) / (6 * h ** 2)
    assert np.allclose(ev, np.sort(ex.ravel()), atol=1e-10)
    assert np.allclose(D, D.T)


# ---------------------------------------------------------------------------
# 2. 热传导: 圆盘中心高斯热源, 与解析解比较, 并检验二阶收敛
# ---------------------------------------------------------------------------
def _disk_poisson_error(h, R=6.0, w=1.0, U=5.0):
    M = NonlocalModel(disk_sample(R, h))
    I = U / (np.pi * w ** 2) * np.exp(-(M.xo ** 2 + M.yo ** 2) / w ** 2)
    n = M.potential(M.thermal_index(I))
    r = np.hypot(M.xo, M.yo)
    exact = U / (2 * np.pi) * (np.log(R / r) - 0.5 * exp1(r ** 2 / w ** 2) + 0.5 * exp1(R ** 2 / w ** 2))
    return np.max(np.abs(n - exact))


def test_poisson_disk_second_order():
    e1, e2 = _disk_poisson_error(0.2), _disk_poisson_error(0.1)
    assert e2 < 2e-3
    assert e1 / e2 > 3.0          # 二阶收敛: 误差比约为 4


# ---------------------------------------------------------------------------
# 3. 对称性传递: 正多边形中折射率的 cos(n theta) 分量 = U/(pi n (n+1)) (r/rho)^n
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n", [3, 4])
def test_symmetry_transfer_coefficient(n):
    rho, U = 10.0, 10.0
    M = NonlocalModel(polygon_sample(n, 0.25, rho=rho))
    I = np.exp(-(M.xo ** 2 + M.yo ** 2) / 2.0)
    I *= U / M.power(np.sqrt(I))
    ng = M.index_grid(M.thermal_index(I))
    r = 3.0
    a, _ = dg.harmonic_amplitudes(ng, M.grid, (0, 0), r, kmax=n)
    assert abs(a[n] / theory.predicted_cos_harmonic(n, U, rho, r) - 1) < 5e-3


# ---------------------------------------------------------------------------
# 4. 牛顿迭代: 基模二次收敛; 定功率模式
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def square_model():
    return NonlocalModel(rectangle_sample(12, 12, 0.3))


def test_newton_fundamental(square_model):
    M = square_model
    w0, _ = guesses.best_guess(M, lambda W: guesses.lg_vortex(M, 0, W), 2.0, np.linspace(1, 3, 5))
    sol = newton.solve(M, w0.real, 2.0)
    assert sol["converged"] and sol["iterations"] <= 10
    h = sol["res_hist"]
    assert h[-1] < 1e-9
    sol_p = newton.solve(M, sol["w"], 1.5, fixed_power=1.05 * sol["power"])
    assert sol_p["converged"]
    assert abs(sol_p["power"] / (1.05 * sol["power"]) - 1) < 1e-8
    assert sol_p["b"] > sol["b"]          # 能流随 b 单调增


def test_vortex_charge_and_irreps():
    M = NonlocalModel(disk_sample(8, 0.25))
    w0, _ = guesses.best_guess(M, lambda W: guesses.lg_vortex(M, 1, W), 2.0, np.linspace(1, 3, 5))
    sol = newton.solve(M, w0, 2.0)
    assert sol["converged"]
    row = dg.mode_summary(M, sol, nsym=4, m=1)
    assert row["total_charge"] == pytest.approx(1.0, abs=1e-3)
    assert row["W_small"] == pytest.approx(1.0, abs=1e-3)
    assert row["mu_dominant"] == 1 and row["mu_weight"] > 0.9999
    assert row["n_satellites"] == 0


# ---------------------------------------------------------------------------
# 5. 稳定性: 基模稳定, 谱的哈密顿对称性; 稠密法与 Arnoldi 法一致
# ---------------------------------------------------------------------------
def test_stability_fundamental(square_model):
    M = square_model
    w0, _ = guesses.best_guess(M, lambda W: guesses.lg_vortex(M, 0, W), 2.0, np.linspace(1, 3, 5))
    sol = newton.solve(M, w0.real, 2.0)
    arn = stability.stability_arnoldi(M, sol, shifts=(0.02,), k=8)
    assert arn["max_growth"] < 1e-6
    dense = stability.stability_dense(M, sol)            # 小网格, 不截断
    lam = dense["eigenvalues"]
    assert dense["max_growth"] < 1e-6
    # Arnoldi 找到的本征值必须出现在稠密谱中
    for l in arn["eigenvalues"]:
        assert np.min(np.abs(lam - l)) < 1e-6


# ---------------------------------------------------------------------------
# 6. 传播: 能流守恒到机器精度, 定态解保持定态
# ---------------------------------------------------------------------------
def test_propagation_conservation(square_model):
    M = square_model
    w0, _ = guesses.best_guess(M, lambda W: guesses.lg_vortex(M, 0, W), 2.0, np.linspace(1, 3, 5))
    sol = newton.solve(M, w0.real, 2.0)
    res = propagate.propagate(M, sol["w"], 5.0, dz=0.02, w_ref=sol["w"], verbose=False)
    assert np.ptp(res["power"]) / res["power"][0] < 1e-10
    assert res["deviation"][-1] < 1e-3


# ---------------------------------------------------------------------------
# 7. 理论选择定则
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("m,n,central,nsat", [
    (2, 3, -1, 3), (3, 4, -1, 4), (2, 2, 0, 2), (3, 2, 1, 2), (4, 2, 0, 4),
    (5, 3, -1, 6), (1, 3, 1, 0), (2, 4, 2, 0), (3, 3, 0, 3), (4, 6, -2, 6)])
def test_selection_rule(m, n, central, nsat):
    r = theory.selection_rule(m, n)
    assert r["central"] == central and r["n_satellites"] == nsat
    assert (r["central"] + sum(k * q for k, q, _ in r["rings"])) == m      # 总荷守恒


def test_selection_rule_rings():
    r = theory.selection_rule(5, 3)
    assert [(k, s) for k, _, s in r["rings"]] == [(3, pytest.approx(1 / 3)), (3, pytest.approx(1.0))]
    assert theory.selection_rule(2, 2)["degenerate"] and not theory.selection_rule(2, 3)["degenerate"]


# ---------------------------------------------------------------------------
# 8. 诊断: 奇点检测与 C_n 投影 (人工构造的场)
# ---------------------------------------------------------------------------
def test_singularity_detection_synthetic():
    M = NonlocalModel(disk_sample(6, 0.1))
    z = M.xo + 1j * M.yo
    env = np.exp(-np.abs(z) ** 2 / 4)
    # 中心 -1 荷 + 三个 +1 卫星 (半径 0.5), 总荷 +2:  (z^3 - a^3) conj(z) / |z|... 用多项式构造
    a = 0.5
    field = (z ** 3 - a ** 3) * np.conj(z) * env          # 零点: z=0 (荷 -1) 与 z^3 = a^3 (三个 +1)
    q = dg.field_grid(M, field)
    sings = dg.find_singularities(q, M.grid, dg.beam_region(q, M.grid) & M.opt_mask)
    st = dg.singularity_structure(sings, (0, 0), 0.05)
    assert st["central_charge"] == -1 and st["n_satellites"] == 3
    assert abs(st["satellite_r_mean"] - a) < 0.01
    w = dg.cn_weights(q, M.grid, 3)
    assert w[-1] > 0.999                                   # m=2 ≡ -1 (mod 3)


def test_multicharge_center_not_aliased():
    """正则 m=3 涡旋: 方格法会混叠成“-1 + 4 卫星”, mode_summary 必须用圆周绕数给出中心荷 +3、无卫星。"""
    M = NonlocalModel(disk_sample(6, 0.2))
    z = M.xo + 1j * M.yo
    sol = dict(w=z ** 3 * np.exp(-np.abs(z) ** 2 / 4), b=1.0, power=1.0)
    row = dg.mode_summary(M, sol)
    assert row["central_charge"] == 3 and row["n_satellites"] == 0 and row["unresolved_charge"] == 0
    assert row["total_charge"] == pytest.approx(3.0, abs=1e-3)


# ---------------------------------------------------------------------------
# 9. 边界温度图样与界面样品
# ---------------------------------------------------------------------------
def test_boundary_temperature_harmonic():
    s = disk_sample(8, 0.2)
    f = harmonic_multipole(3, 0.5, 8)
    s.set_boundary_temperature(f)
    M = NonlocalModel(s)
    n_bc = M.potential(M.n_bc)
    assert np.max(np.abs(n_bc - f(M.xo, M.yo))) < 2e-3


def test_surface_sample_soliton():
    M = NonlocalModel(surface_interface_sample(L=16, nd=-100, h=0.4))
    w0 = guesses.best_guess(M, lambda W: guesses.multipole_ring(M, 1, 0.0, W, center=(-1.5, 0)),
                            2.0, np.linspace(1, 3, 5))[0]
    sol = newton.solve(M, w0.real, 2.0)
    assert sol["converged"]
    xc, _ = dg.center_of_mass(M, sol["w"])
    assert -4 < xc < 0                                     # 光束贴在热介质一侧的界面附近
