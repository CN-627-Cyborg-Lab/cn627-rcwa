"""
1D-RCWA engine for binary gratings (TE, normal incidence).
一维RCWA引擎，适用于二元光栅（TE偏振，正入射）。

Algorithm / 算法:
    Moharam et al., JOSA A 12(5), 1068 (1995).
S-matrix cascade / S矩阵级联:
    Li, JOSA A 13(5), 1024 (1996).

Limitations / 限制:
    - TE polarization only / 仅TE偏振（TM需要逆规则）
    - Normal incidence only / 仅正入射（θ=0）
    - Single grating layer / 单层光栅
    - Real or complex permittivity / 实数或复数介电常数
"""

import numpy as np
from .materials import get_si_eps


# ── Branch-cut fix / 分支切割修正 ──────────────────────
def _fix_gamma(g):
    """
    Ensure propagating modes have Re(γ)>0, evanescent have Im(γ)>0.
    Regularize near-zero values to avoid singular matrices at grazing orders.
    确保传播模式 Re(γ)>0，倏逝模式 Im(γ)>0。
    正则化近零值，避免掠射级次导致奇异矩阵。
    """
    g = np.asarray(g, dtype=np.complex128)
    g = np.where(np.imag(g) < 0, -g, g)
    g = np.where(
        (np.abs(np.imag(g)) < 1e-10) & (np.real(g) < 0),
        -g, g
    )
    # Regularize: when λ/Λ is integer, some orders have γ exactly 0.
    # These are grazing modes carrying no power. A tiny offset avoids
    # singular matrices without affecting physical results.
    # 正则化：当 λ/Λ 为整数时，某些级次的 γ 恰好为0。
    # 这些是不携带能量的掠射模式。微小偏移避免奇异矩阵，不影响物理结果。
    g = np.where(np.abs(g) < 1e-12, 1e-12 + 0j, g)
    return g


# ── Redheffer star product / Redheffer星积（Li1996 Eq.7）──
def _redheffer(SA, SB):
    """
    Cascade two S-matrices via Redheffer star product.
    通过Redheffer星积级联两个S矩阵。
    """
    SA11, SA12, SA21, SA22 = SA
    SB11, SB12, SB21, SB22 = SB
    N = SA11.shape[0]
    I = np.eye(N, dtype=np.complex128)

    F1 = np.linalg.inv(I - SB11 @ SA22)
    F2 = np.linalg.inv(I - SA22 @ SB11)

    C11 = SA11 + SA12 @ F1 @ SB11 @ SA21
    C12 = SA12 @ F1 @ SB12
    C21 = SB21 @ F2 @ SA21
    C22 = SB22 + SB21 @ F2 @ SA22 @ SB12
    return C11, C12, C21, C22


def _S_identity(N):
    """
    Identity S-matrix (transparent layer).
    单位S矩阵（透明层）。
    """
    I = np.eye(N, dtype=np.complex128)
    Z = np.zeros((N, N), dtype=np.complex128)
    return Z.copy(), I.copy(), I.copy(), Z.copy()


# ── Toeplitz matrix / Toeplitz矩阵（Moharam1995 Eq.3-5）──
def _build_toeplitz(ff, eps_high, eps_low, n_orders):
    """
    Build permittivity Fourier coefficient Toeplitz matrix.
    构建介电常数傅里叶系数的Toeplitz矩阵。
    """
    N = 2 * n_orders + 1
    center = 2 * n_orders
    m = np.arange(-2 * n_orders, 2 * n_orders + 1, dtype=float)

    eps_hat = np.zeros(4 * n_orders + 1, dtype=np.complex128)
    eps_hat[center] = ff * eps_high + (1.0 - ff) * eps_low

    nz = m != 0
    eps_hat[nz] = (
        (eps_high - eps_low)
        * np.sin(m[nz] * np.pi * ff) / (m[nz] * np.pi)
    )

    idx = np.arange(N)
    diff = idx[:, None] - idx[None, :]
    return eps_hat[diff + center]


# ── Interface S-matrix / 界面S矩阵（Moharam1995 Eq.17-20）──
def _interface_S(gamma_L, W_L, gamma_R, W_R):
    """
    S-matrix at interface between two regions (TE).
    两个区域界面处的S矩阵（TE偏振）。
    """
    N = len(gamma_L)
    I = np.eye(N, dtype=np.complex128)

    F = W_L @ np.diag(gamma_L)
    G = W_R @ np.diag(gamma_R)

    WR_inv = np.linalg.inv(W_R)
    WL_inv = np.linalg.inv(W_L)

    M = np.linalg.inv(F) @ G @ WR_inv @ W_L
    S11 = np.linalg.inv(I + M) @ (I - M)
    S21 = WR_inv @ W_L @ (I + S11)

    M_rev = np.linalg.inv(G) @ F @ WL_inv @ W_R
    S22 = np.linalg.inv(I + M_rev) @ (I - M_rev)
    S12 = WL_inv @ W_R @ (I + S22)

    return S11, S12, S21, S22


# ── Propagation S-matrix / 传播S矩阵（Moharam1995 Eq.22-23）──
def _propagation_S(gamma, k0, d):
    """
    S-matrix for propagation through a uniform layer.
    均匀层内传播的S矩阵。
    """
    X = np.exp(1j * gamma * k0 * d)
    N = len(gamma)
    Z = np.zeros((N, N), dtype=np.complex128)
    return Z.copy(), np.diag(X), np.diag(X), Z.copy()


# ── Full grating layer S-matrix / 完整光栅层S矩阵 ──────
def _grating_S(kx_norm, eps_toep, k0, d, eps_inc, eps_sub, ff, eps_high, eps_low):
    """
    Build complete S-matrix: incidence → grating → substrate.
    构建完整S矩阵：入射层 → 光栅层 → 基底层。
    """
    N = eps_toep.shape[0]

    # --- 新增的均匀层特殊处理（跳过特征值求解，防止奇异矩阵） ---
    if ff == 1.0 or eps_high == eps_low:
        eps_uniform = eps_high
    elif ff == 0.0:
        eps_uniform = eps_low
    else:
        eps_uniform = None

    if eps_uniform is not None:
        # 均匀介质直接使用解析解，特征矩阵 W 退化为单位阵
        gamma_g = _fix_gamma(
            np.sqrt((complex(eps_uniform) - kx_norm ** 2).astype(np.complex128))
        )
        W = np.eye(N, dtype=np.complex128)
    else:
        # 周期性光栅，执行常规数值求解
        eigvals, W = np.linalg.eig(eps_toep - np.diag(kx_norm ** 2))
        gamma_g = _fix_gamma(np.sqrt(eigvals.astype(np.complex128)))
    # -----------------------------------------------------------

    gamma_inc = _fix_gamma(
        np.sqrt((complex(eps_inc) - kx_norm ** 2).astype(np.complex128))
    )
    gamma_sub = _fix_gamma(
        np.sqrt((complex(eps_sub) - kx_norm ** 2).astype(np.complex128))
    )

    I_mat = np.eye(N, dtype=np.complex128)

    S_in = _interface_S(gamma_inc, I_mat, gamma_g, W)
    S_pr = _propagation_S(gamma_g, k0, d)
    S_ou = _interface_S(gamma_g, W, gamma_sub, I_mat)

    S = _S_identity(N)
    S = _redheffer(S, S_in)
    S = _redheffer(S, S_pr)
    S = _redheffer(S, S_ou)
    return S, gamma_inc, gamma_sub


# ── Public API / 公开接口 ──────────────────────────────
def rcwa_te(wl_nm, period_nm, height_nm, ff,
            eps_high=None, eps_low=1.0,
            eps_inc=1.0, eps_sub=1.0,
            n_orders=10):
    """
    1D-RCWA for a single binary grating layer (TE, normal incidence).
    单层二元光栅的一维RCWA计算（TE偏振，正入射）。

    Parameters / 参数
    ----------
    wl_nm : float
        Wavelength in nm. / 波长（纳米）。
    period_nm : float
        Grating period in nm. / 光栅周期（纳米）。
    height_nm : float
        Grating layer thickness in nm. / 光栅层厚度（纳米）。
    ff : float
        Fill fraction of high-index material. / 高折射率材料填充比。
    eps_high : complex or None
        Permittivity of high-index region. None = Si dispersion.
        高折射率区域介电常数。None 则使用硅色散数据。
    eps_low : complex
        Permittivity of low-index region (default: Air).
        低折射率区域介电常数（默认：空气）。
    eps_inc : complex
        Incidence medium permittivity (default: Air).
        入射介质介电常数（默认：空气）。
    eps_sub : complex
        Substrate permittivity (default: Air).
        基底介电常数（默认：空气）。
    n_orders : int
        Truncation order. Total modes = 2*n_orders + 1.
        截断阶数。总模式数 = 2*n_orders + 1。

    Returns / 返回
    -------
    R : float
        Total reflectance. / 总反射率。
    T : float
        Total transmittance. / 总透射率。
    A : float
        Absorptance = max(0, 1-R-T). / 吸收率。
    energy_err : float
        Energy conservation error |R+T-1| (lossless) or |R+T+A-1|.
        能量守恒误差。
    """
    lam = wl_nm * 1e-9
    d = height_nm * 1e-9
    Lam = period_nm * 1e-9
    k0 = 2.0 * np.pi / lam
    N = 2 * n_orders + 1

    m = np.arange(-n_orders, n_orders + 1, dtype=float)
    kx_norm = m * lam / Lam

    if eps_high is None:
        eps_high = get_si_eps(wl_nm)

    E = _build_toeplitz(ff, eps_high, eps_low, n_orders)
    S, gamma_inc, gamma_sub = _grating_S(kx_norm, E, k0, d, eps_inc, eps_sub, ff, eps_high, eps_low)
    # Incident field: zeroth order only
    # 入射场：仅零级
    e_inc = np.zeros(N, dtype=np.complex128)
    e_inc[n_orders] = 1.0

    r = S[0] @ e_inc    # S11 @ e_inc → reflected amplitudes / 反射振幅
    t = S[2] @ e_inc    # S21 @ e_inc → transmitted amplitudes / 透射振幅

    g0 = np.real(gamma_inc[n_orders])
    prop_inc = np.real(gamma_inc) > 1e-10
    prop_sub = np.real(gamma_sub) > 1e-10

    R = np.sum(np.abs(r[prop_inc]) ** 2 * np.real(gamma_inc[prop_inc])) / g0
    T = np.sum(np.abs(t[prop_sub]) ** 2 * np.real(gamma_sub[prop_sub])) / g0

    energy_err = abs(R + T - 1.0)
    A = max(0.0, 1.0 - R - T)

    return R, T, A, energy_err


def rcwa_spectrum(wl_array, period_nm, height_nm, ff, **kwargs):
    """
    Run rcwa_te over an array of wavelengths.
    对波长数组批量运行 rcwa_te。

    Returns dict with keys: 'wl', 'R', 'T', 'A', 'energy_err'.
    返回字典，键为：'wl', 'R', 'T', 'A', 'energy_err'。
    """
    results = {k: [] for k in ('wl', 'R', 'T', 'A', 'energy_err')}
    for wl in wl_array:
        R, T, A, err = rcwa_te(float(wl), period_nm, height_nm, ff, **kwargs)
        results['wl'].append(float(wl))
        results['R'].append(R)
        results['T'].append(T)
        results['A'].append(A)
        results['energy_err'].append(err)

    return {k: np.array(v) for k, v in results.items()}