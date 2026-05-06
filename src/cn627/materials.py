"""
Material database for CN-627 RCWA.
CN-627 RCWA 材料数据库。

Si dispersion: M. A. Green, Solar Energy Mat. & Solar Cells 92(11),
1305-1310 (2008). Range: 250-1200 nm.
硅色散数据：Green 2008，波长范围 250-1200 nm。

Ge/Se/SiO2: constant n values near 10.6 µm from Niraula et al.,
Opt. Express 22(21), 25817 (2014).
Ge/Se/SiO2：10.6 µm 附近常数折射率，来自 Niraula 2014。
"""

import numpy as np

# --- Si dispersion data (Green 2008) ---
# --- 硅色散数据（Green 2008）---
_SI_WL = np.array([
    250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
    750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200
], dtype=float)

_SI_N = np.array([
    1.665, 4.976, 5.494, 5.613, 4.691, 4.294, 4.077,
    3.940, 3.844, 3.772, 3.717, 3.675, 3.641, 3.614,
    3.591, 3.572, 3.556, 3.542, 3.530, 3.520
], dtype=float)

_SI_K = np.array([
    3.665, 4.234, 2.938, 0.296, 0.086302, 0.044165,
    0.027968, 0.019934, 0.014431, 0.010528, 0.0078185,
    0.0054113, 0.0036120, 0.0021701, 0.0011793,
    0.0005093, 0.0001362, 3.0637e-05, 6.2230e-06,
    2.1008e-07
], dtype=float)


def get_si_nk(wl_nm):
    """
    Get Si refractive index (n, k) at given wavelength (nm).
    获取指定波长（nm）处硅的折射率 (n, k)。
    """
    n = float(np.interp(wl_nm, _SI_WL, _SI_N))
    k = float(np.interp(wl_nm, _SI_WL, _SI_K))
    return n, k


def get_si_eps(wl_nm):
    """
    Get Si complex permittivity at given wavelength (nm).
    获取指定波长（nm）处硅的复介电常数。
    """
    n, k = get_si_nk(wl_nm)
    return complex(n, k) ** 2


# --- Constant-index materials (mid-IR) ---
# --- 常数折射率材料（中红外）---
MATERIALS = {
    "Ge":   {"n": 4.00, "eps": 16.0000},
    "Se":   {"n": 2.64, "eps":  6.9696},
    "SiO2": {"n": 1.40, "eps":  1.9600},
    "Air":  {"n": 1.00, "eps":  1.0000},
}


def get_eps_const(name):
    """
    Get constant permittivity by material name.
    通过材料名称获取常数介电常数。
    """
    return MATERIALS[name]["eps"]