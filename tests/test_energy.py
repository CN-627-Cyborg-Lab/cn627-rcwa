"""
Energy conservation tests.
能量守恒测试。

Key physics / 关键物理:
  - Lossless materials: R + T = 1 (strict)
    无损耗材料：R + T = 1（严格）
  - Absorbing materials: R + T < 1, R + T + A = 1
    吸收材料：R + T < 1，R + T + A = 1
  - In ALL cases: R + T must NEVER exceed 1 (no energy creation)
    所有情况：R + T 绝不能超过1（不能凭空产生能量）
"""
import numpy as np
import pytest
from cn627 import rcwa_te
from cn627.materials import get_si_nk, MATERIALS


def test_si_grating_energy():
    """
    Si/Air grating: R+T must never exceed 1.
    For absorbing wavelengths (400-500nm), R+T < 1 is correct physics.
    For near-lossless wavelengths (>900nm), R+T ≈ 1.
    硅/空气光栅：R+T绝不能超过1。
    吸收波段（400-500nm）R+T<1是正确物理。
    近无损波段（>900nm）R+T≈1。
    """
    for wl in np.arange(400, 1001, 50):
        R, T, A, err = rcwa_te(
            wl, period_nm=700.0, height_nm=200.0, ff=0.5,
            n_orders=15
        )
        # No energy creation (universal)
        # 不能凭空产生能量（普适）
        assert R + T <= 1.0 + 1e-6, f"R+T > 1 at {wl} nm (non-physical)"

        # For near-lossless wavelengths, check conservation
        # 对近无损波段，检查守恒
        _, k = get_si_nk(wl)
        if k < 1e-3:
            assert err < 0.01, f"Energy error {err:.2e} at {wl} nm (low-loss)"


def test_lossless_grating_energy():
    """
    Lossless Ge/Se grating: R+T must equal 1.
    无损耗Ge/Se光栅：R+T必须等于1。
    """
    eps_ge = MATERIALS["Ge"]["eps"]
    eps_se = MATERIALS["Se"]["eps"]
    eps_sub = MATERIALS["SiO2"]["eps"]

    for wl in np.arange(9000, 12001, 200):
        R, T, A, err = rcwa_te(
            wl, period_nm=6910.0, height_nm=3800.0, ff=0.42,
            eps_high=eps_ge, eps_low=eps_se,
            eps_sub=eps_sub, n_orders=21
        )
        assert err < 1e-6, f"|R+T-1| = {err:.2e} at {wl} nm"