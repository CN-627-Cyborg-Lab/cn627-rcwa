"""
Example: Si/Air binary grating spectral scan.
示例：硅/空气二元光栅光谱扫描。

Demonstrates basic usage of cn627 for a visible-range Si grating.
展示cn627在可见光范围硅光栅上的基本用法。
"""
import numpy as np
from cn627 import rcwa_spectrum

# --- Parameters / 参数 ---
wl = np.arange(400, 1001, 5)   # nm / 纳米
period = 700.0    # nm / 纳米
height = 200.0    # nm / 纳米
ff = 0.5

# --- Run / 运行 ---
result = rcwa_spectrum(wl, period, height, ff, n_orders=21)

# --- Print summary / 输出摘要 ---
i_peak = np.argmax(result['T'])
print(f"T peak / T峰值: {result['T'][i_peak]:.4f} @ {result['wl'][i_peak]:.0f} nm")
print(f"Max energy error / 最大能量误差: {result['energy_err'].max():.2e}")
print(f"Points with R+T>1 / R+T>1的点数: "
      f"{np.sum(result['R'] + result['T'] > 1 + 1e-6)}")

# --- Save CSV / 保存CSV ---
np.savetxt("si_grating_spectrum.csv",
           np.column_stack([result['wl'], result['R'],
                            result['T'], result['A']]),
           delimiter=",", header="wavelength_nm,R,T,A",
           comments="")
print("Saved / 已保存: si_grating_spectrum.csv")