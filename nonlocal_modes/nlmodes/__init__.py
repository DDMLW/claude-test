"""
nlmodes —— 热非局域非线性介质中特殊模态(涡旋/多极/表面态)的数值工具包

主要模块:
    geometry    样品几何与网格
    operators   稀疏差分算子
    model       非局域模型 (折射率 = 热传导方程的解)
    newton      定态解牛顿迭代与参数延拓
    stability   线性稳定性分析
    propagate   传播动力学 (Crank-Nicolson 分步法)
    diagnostics 相位奇点、C_n 不可约表示等诊断
    guesses     初始猜测
    theory      理论公式 (选择定则、各向异性系数等)
    io_utils    结果保存/读取
    plotting    标准作图
"""

__version__ = "0.1.0"
