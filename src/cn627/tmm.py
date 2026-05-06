"""
Transfer Matrix Method — single-slab analytical solution.
传递矩阵法 - 单层薄膜解析解。

Used as ground-truth reference for RCWA flat-slab validation (FF=1.0).
用作RCWA平板验证（FF=1.0）的基准参考。

Supports complex permittivity (absorbing films).
支持复数介电常数（吸收薄膜）。
"""

import numpy as np


def tmm_slab(wl_nm, thickness_nm, eps_slab,
             eps_inc=1.0, eps_sub=1.0):
    """
    TMM for a single slab (TE, normal incidence).
    单层薄膜的TMM计算（TE偏振，正入射）。

    Parameters / 参数
    ----------
    wl_nm : float
        Wavelength in nm. / 波长（纳米）。
    thickness_nm : float
        Slab thickness in nm. / 薄膜厚度（纳米）。
    eps_slab : complex
        Slab permittivity. / 薄膜介电常数。
    eps_inc : float
        Incidence medium permittivity. / 入射介质介电常数。
    eps_sub : float
        Substrate permittivity. / 基底介电常数。

    Returns / 返回
    -------
    R, T, A : float
        Reflectance, transmittance, absorptance.
        反射率、透射率、吸收率。
    """
    lam = wl_nm * 1e-9
    d = thickness_nm * 1e-9
    k0 = 2.0 * np.pi / lam

    n_inc = np.sqrt(complex(eps_inc))
    n1 = np.sqrt(complex(eps_slab))
    n_sub = np.sqrt(complex(eps_sub))

    # Branch-cut fix / 分支切割修正
    for n in (n1, n_sub):
        if np.real(n) < 0:
            n = -n  # noqa

    kz0 = n_inc * k0
    kz1 = n1 * k0
    kz2 = n_sub * k0

    # Fresnel coefficients / 菲涅尔系数
    r01 = (kz0 - kz1) / (kz0 + kz1)
    t01 = 2 * kz0 / (kz0 + kz1)
    r12 = (kz1 - kz2) / (kz1 + kz2)
    t12 = 2 * kz1 / (kz1 + kz2)

    # Airy formula / 艾里公式
    phi = kz1 * d
    ep = np.exp(2j * phi)
    den = 1 + r01 * r12 * ep

    r = (r01 + r12 * ep) / den
    t = t01 * t12 * np.exp(1j * phi) / den

    R = abs(r) ** 2
    T = abs(t) ** 2 * np.real(kz2) / np.real(kz0)
    A = max(0.0, 1.0 - R - T)
    return R, T, A