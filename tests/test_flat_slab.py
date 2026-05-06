"""
Flat slab validation: RCWA(FF=1) must match TMM.
平板验证：RCWA（FF=1）必须与TMM匹配。
"""
import numpy as np
import pytest
from cn627 import rcwa_te, tmm_slab
from cn627.materials import get_si_eps, MATERIALS


@pytest.mark.parametrize("wl", [400, 500, 600, 700, 800, 900, 1000])
def test_si_slab_vs_tmm(wl):
    """
    Si slab in air, 200 nm thick.
    空气中200nm硅薄膜。
    """
    eps = get_si_eps(wl)
    R_tmm, T_tmm, _ = tmm_slab(wl, 200.0, eps)
    R_rcwa, T_rcwa, _, err = rcwa_te(
        wl, period_nm=1000.0, height_nm=200.0, ff=1.0,
        eps_high=eps, eps_low=eps, n_orders=10
    )
    assert abs(R_rcwa - R_tmm) / max(R_tmm, 1e-12) < 0.01
    assert abs(T_rcwa - T_tmm) / max(T_tmm, 1e-12) < 0.01


@pytest.mark.parametrize("wl", [10000, 10400, 10600, 10800, 11200])
def test_ge_slab_on_sio2(wl):
    """
    Ge slab on SiO2 substrate, 3800 nm thick.
    SiO2基底上3800nm锗薄膜。
    """
    eps_ge = MATERIALS["Ge"]["eps"]
    eps_sub = MATERIALS["SiO2"]["eps"]
    R_tmm, T_tmm, _ = tmm_slab(wl, 3800.0, eps_ge, eps_sub=eps_sub)
    R_rcwa, T_rcwa, _, err = rcwa_te(
        wl, period_nm=9000.0, height_nm=3800.0, ff=1.0,
        eps_high=eps_ge, eps_low=eps_ge, eps_sub=eps_sub, n_orders=10
    )
    assert abs(R_rcwa - R_tmm) / max(R_tmm, 1e-12) < 0.01
    assert abs(T_rcwa - T_tmm) / max(T_tmm, 1e-12) < 0.01