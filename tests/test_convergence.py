"""
Convergence test: results must stabilize with increasing orders.
收敛测试：结果必须随截断阶数增加而稳定。

Note: period chosen to avoid λ/Λ = integer, which causes
grazing-order singularity.
注意：周期选择避开 λ/Λ 为整数，防止掠射级次奇异。
"""
import numpy as np
from cn627 import rcwa_te


def test_convergence_si():
    """
    R, T should converge as n_orders increases.
    Avoid λ/Λ = integer to prevent grazing-order singularity.
    Check convergence only at higher orders (n >= 31) where
    results are expected to stabilize.
    R、T应随n_orders增加而收敛。
    避开 λ/Λ 为整数，防止掠射级次奇异。
    仅在高阶（n>=31）检查收敛，低阶尚未稳定属正常现象。
    """
    prev_R, prev_T = None, None
    for n in [5, 11, 21, 31, 41, 51]:
        R, T, _, _ = rcwa_te(
            700.0, period_nm=600.0, height_nm=200.0, ff=0.5,
            n_orders=n
        )
        if prev_R is not None and n >= 31:
            assert abs(R - prev_R) < 1e-3, f"R not converged at n={n}"
            assert abs(T - prev_T) < 1e-3, f"T not converged at n={n}"
        prev_R, prev_T = R, T