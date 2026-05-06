"""
Validation: reproduce Device 1 from Niraula et al., Opt. Express 22, 25817 (2014).
验证：复现 Niraula 2014 论文中的 Device 1。

Structure / 结构:
    Ge/Se binary grating (Λ=6910nm, d=3800nm, FF=0.42) on SiO2.
    SiO2基底上的Ge/Se二元光栅（Λ=6910nm, d=3800nm, FF=0.42）。
Expected / 预期:
    Narrow transmission peak near 10.6 µm.
    10.6 µm 附近出现窄透射峰。
"""
import numpy as np
from cn627 import rcwa_te, rcwa_spectrum
from cn627.materials import MATERIALS

# --- Device 1 parameters / Device 1 参数 ---
PERIOD = 6910.0   # nm / 纳米
HEIGHT = 3800.0   # nm / 纳米
FF = 0.42
EPS_GE = MATERIALS["Ge"]["eps"]      # 16.0
EPS_SE = MATERIALS["Se"]["eps"]      # 6.9696
EPS_SUB = MATERIALS["SiO2"]["eps"]   # 1.96

common = dict(
    period_nm=PERIOD, height_nm=HEIGHT, ff=FF,
    eps_high=EPS_GE, eps_low=EPS_SE, eps_sub=EPS_SUB
)

# --- Step 1: Convergence / 第一步：收敛性测试 ---
print("=== Convergence test @ 10600 nm / 收敛测试 @ 10600 nm ===")
prev = None
for n in [5, 11, 21, 31, 41, 55]:
    R, T, A, err = rcwa_te(10600.0, n_orders=n, **common)
    delta = (f"ΔR={abs(R-prev[0]):.1e} ΔT={abs(T-prev[1]):.1e}"
             if prev else "—")
    print(f"  n={n:3d}  R={R:.6f}  T={T:.6f}  err={err:.1e}  {delta}")
    prev = (R, T)

# --- Step 2: Broadband scan / 第二步：宽带扫描 ---
print("\n=== Broadband scan 8-14 µm / 宽带扫描 8-14 µm ===")
wl_broad = np.arange(8000, 14001, 10)
res = rcwa_spectrum(wl_broad, n_orders=41, **common)

i_peak = np.argmax(res['T'])
print(f"T peak / T峰值: {res['T'][i_peak]:.4f} @ {res['wl'][i_peak]:.0f} nm")
print(f"Max energy error / 最大能量误差: {res['energy_err'].max():.2e}")

# --- Step 3: Fine scan around peak / 第三步：峰值附近精细扫描 ---
print("\n=== Fine scan / 精细扫描 ===")
center = res['wl'][i_peak]
wl_fine = np.arange(center - 50, center + 50, 0.1)
res_f = rcwa_spectrum(wl_fine, n_orders=41, **common)

i_fine = np.argmax(res_f['T'])
peak_wl = res_f['wl'][i_fine]
peak_T = res_f['T'][i_fine]
print(f"Fine peak / 精细峰值: T={peak_T:.6f} @ {peak_wl:.1f} nm")

# --- Step 4: Paper comparison / 第四步：论文对比 ---
print("\n=== Comparison with Niraula 2014 / 与Niraula 2014对比 ===")
print(f"  Paper / 论文:   λ ≈ 10600 nm,  T ≈ 0.90")
print(f"  Ours / 本项目:  λ = {peak_wl:.1f} nm,  T = {peak_T:.4f}")
print(f"  Δλ = {peak_wl - 10600:+.1f} nm,  ΔT = {peak_T - 0.90:+.4f}")
print(f"\n  Note: deviations expected due to constant-n model (no dispersion).")
print(f"  注意：偏差来自常数折射率模型（未考虑色散）。")