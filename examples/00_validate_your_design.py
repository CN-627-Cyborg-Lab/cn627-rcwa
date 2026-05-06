"""
=============================================================================
CN-627 RCWA Validator — Design Validation Script
CN-627 RCWA 验证器 — 设计验证脚本

HOW TO USE / 使用方法:
    1. Edit Section 1: fill in your design parameters
       修改第 1 节：填入你的设计参数
    2. Edit Section 2: choose wavelength mode
       修改第 2 节：选择波长模式
    3. Run:  python examples/00_validate_your_design.py
       运行：python examples/00_validate_your_design.py
    4. Read the report in your terminal
       在终端阅读报告

STRUCTURE DIAGRAM / 结构示意图:

    ┌─────────────────────────────┐
    │  Incidence medium (eps_inc) │  ← light enters here / 光从这里入射
    │  default: Air, n=1.0        │
    ├──────────┬─────┬─────┬──────┤
    │ eps_high │ low │ high│ low  │  ← grating layer / 光栅层
    │          │     │     │      │     period Λ, height H, fill fraction FF
    ├──────────┴─────┴─────┴──────┤
    │  Substrate (eps_sub)        │  ← light exits here / 光从这里透射
    │  default: Air, n=1.0        │
    └─────────────────────────────┘

MATERIAL CONVERSION / 材料转换:
    eps = complex(n, k) ** 2
    e.g. n=3.48, k=0.00548 → eps = (3.48+0.00548j)² = 12.11+0.038j
    If no absorption / 无吸收：k=0, eps = n**2
    e.g. n=1.47 → eps = 1.47**2 = 2.1609

=============================================================================
"""
import sys
import os
import time
import numpy as np

# Add project root to path / 添加项目根目录到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from cn627 import rcwa_te, rcwa_spectrum
from cn627.materials import get_si_eps


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: DESIGN PARAMETERS — EDIT HERE                                 ║
# ║  第1节：设计参数 — 在这里修改                                               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

DESIGN_NAME = "My Design"   # give your design a name / 给设计起个名字

# --- Grating geometry / 光栅几何参数 ---
PERIOD = 500.0           # nm, grating period Λ / 光栅周期
HEIGHT = 600.0           # nm, grating height H / 光栅高度
FF = 0.5                 # fill fraction (0 to 1) / 占空比（0到1）

# --- Materials / 材料 ---
#
# How to convert / 如何转换:
#   refractive index n, k  →  eps = complex(n, k) ** 2
#   折射率 n, k            →  介电常数 eps = complex(n, k) ** 2
#
# Examples / 示例:
#   Air    / 空气:    n=1.00          →  eps = 1.0
#   SiO2   / 二氧化硅: n=1.47         →  eps = 1.47**2 = 2.1609
#   Si     / 硅:      n=3.48, k=0.005 →  eps = (3.48+0.005j)**2
#   Ge     / 锗:      n=4.00          →  eps = 16.0
#   Se     / 硒:      n=2.64          →  eps = 6.9696
#
# Special: set EPS_HIGH = None to use built-in Si dispersion (250-1200 nm)
# 特殊：设 EPS_HIGH = None 使用内置硅色散数据（250-1200 nm）

EPS_HIGH = complex(3.48, 0.00548) ** 2   # high-index region / 高折射率区域
EPS_LOW  = 1.0                           # low-index region / 低折射率区域
EPS_INC  = 1.0                           # incidence medium / 入射介质
EPS_SUB  = 1.47 ** 2                     # substrate / 基底


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: WAVELENGTH SELECTION — CHOOSE ONE MODE                        ║
# ║  第2节：波长选择 — 选择一种模式                                             ║
# ╠═══════════════════════════════════════════════════════════════════════════╣
# ║                                                                           ║
# ║  Mode A: specific wavelengths    / 模式A：指定波长点                       ║
# ║  Mode B: range scan              / 模式B：范围扫描                         ║
# ║  Mode C: both (specific + scan)  / 模式C：两者都要                         ║
# ║                                                                           ║
# ║  Set WL_MODE = "A", "B", or "C" below.                                    ║
# ║  在下方设置 WL_MODE = "A", "B", 或 "C"。                                   ║
# ║                                                                           ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

WL_MODE = "C"   # "A" = specific, "B" = scan, "C" = both

# --- Mode A: specific wavelengths / 模式A：指定波长点 ---
# List the exact wavelengths you want to test (nm)
# 列出你想测试的精确波长（nm）
WL_SPECIFIC = [850, 880, 900, 905, 910, 950, 1000]

# --- Mode B: range scan / 模式B：范围扫描 ---
SCAN_START = 900.0    # nm, start / 起点
SCAN_END   = 910.0    # nm, end / 终点
SCAN_STEP  = 0.2      # nm, step / 步长


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: COMPUTATION SETTINGS — OPTIONAL                               ║
# ║  第3节：计算设置 — 可选                                                    ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

N_ORDERS = 45           # truncation order / 截断阶数
                        # 15 = fast / 快速
                        # 45 = accurate (recommended) / 精确（推荐）
                        # 65 = overkill / 过度

RUN_CONVERGENCE = True  # True = run convergence test / 运行收敛测试
                        # False = skip (faster) / 跳过（更快）

SAVE_CSV = True         # True = save results to CSV / 保存CSV
CSV_FILENAME = "validation_result.csv"


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: VALIDATION ENGINE — DO NOT EDIT BELOW                         ║
# ║  第4节：验证引擎 — 以下内容无需修改                                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

def fmt_eps(eps):
    """Format complex eps for display / 格式化复数介电常数"""
    if isinstance(eps, (int, float)):
        return f"{eps:.4f}"
    elif isinstance(eps, complex):
        if abs(eps.imag) < 1e-10:
            return f"{eps.real:.4f}"
        return f"{eps.real:.4f}{eps.imag:+.4f}j"
    return str(eps)


def annotate(R, T, A, cons):
    """Generate annotation for a data row / 生成数据行注释"""
    if R + T > 1.0 + 1e-6:
        return "❌ R+T>1 非物理"
    if cons > 1e-3:
        return "⚠️  守恒误差大"
    if A >= 0.8:
        return "★★ 极高吸收"
    if A >= 0.5:
        return "★ 高吸收"
    if A >= 0.3:
        return "△ 中吸收"
    if A >= 0.1:
        return "○ 低吸收"
    return ""


def print_header(title):
    """Print section header / 打印章节标题"""
    print("\n" + "=" * 72)
    print(f"  {title}")
    print("=" * 72)


def print_table(wl_list, results, mark_peak=True):
    """
    Print formatted results table / 打印格式化结果表格
    results: list of (wl, R, T, A, cons)
    """
    # Find peak A index
    A_vals = [r[3] for r in results]
    i_peak = int(np.argmax(A_vals))

    print(f"\n  {'λ(nm)':>8} │ {'R':>9} {'T':>9} {'A':>9} │"
          f" {'R+T+A':>10} {'守恒误差':>10} │ {'备注'}")
    print("  " + "─" * 72)

    for i, (wl, R, T, A, cons) in enumerate(results):
        note = annotate(R, T, A, cons)
        if mark_peak and i == i_peak:
            note = "◄ PEAK " + note
        print(f"  {wl:>8.1f} │ {R:>9.5f} {T:>9.5f} {A:>9.5f} │"
              f" {R+T+A:>10.7f} {cons:>10.2e} │ {note}")

    print("  " + "─" * 72)
    return i_peak


def run():
    """Main validation routine / 主验证流程"""
    t_start = time.time()

    # ── Title / 标题 ──
    print("\n╔" + "═" * 70 + "╗")
    print("║" + " CN-627 RCWA Validator — Design Validation Report ".center(70) + "║")
    print("║" + " CN-627 RCWA 验证器 — 设计验证报告 ".center(64) + "║")
    print("╚" + "═" * 70 + "╝")

    # ── Step 1: Parameters / 参数确认 ──
    print_header("Step 1: Design Parameters / 设计参数")
    print(f"""
  Design name / 设计名称:    {DESIGN_NAME}
  Period / 周期 Λ:           {PERIOD:.1f} nm
  Height / 高度 H:           {HEIGHT:.1f} nm
  Fill fraction / 占空比 FF: {FF:.4f}  (linewidth / 线宽 = {PERIOD * FF:.1f} nm)
  eps_high:                  {fmt_eps(EPS_HIGH)}{"  (auto Si dispersion)" if EPS_HIGH is None else ""}
  eps_low:                   {fmt_eps(EPS_LOW)}
  eps_inc:                   {fmt_eps(EPS_INC)}
  eps_sub:                   {fmt_eps(EPS_SUB)}
  n_orders:                  {N_ORDERS}  ({2*N_ORDERS+1} modes / 模式)
  Wavelength mode:           {WL_MODE}""")

    # Build common kwargs / 构建公共参数
    common = dict(
        period_nm=PERIOD, height_nm=HEIGHT, ff=FF,
        eps_low=EPS_LOW, eps_inc=EPS_INC, eps_sub=EPS_SUB,
        n_orders=N_ORDERS
    )
    if EPS_HIGH is not None:
        common['eps_high'] = EPS_HIGH

    # ── Step 2: Convergence / 收敛测试 ──
    if RUN_CONVERGENCE:
        # Pick a test wavelength
        if WL_MODE in ("A", "C") and WL_SPECIFIC:
            test_wl = WL_SPECIFIC[len(WL_SPECIFIC) // 2]
        else:
            test_wl = (SCAN_START + SCAN_END) / 2

        print_header(f"Step 2: Convergence Test @ {test_wl:.1f} nm / 收敛测试")

        conv_kw = dict(common)
        conv_kw.pop('n_orders', None)

        print(f"\n  {'n':>5} {'modes':>6} {'R':>10} {'T':>10} {'A':>10}"
              f" {'ΔR':>10} {'ΔT':>10} {'守恒':>10}")
        print("  " + "─" * 68)

        prev_R, prev_T = None, None
        for n in [5, 11, 21, 31, 45, 65]:
            Ri, Ti, Ai, ei = rcwa_te(wl_nm=test_wl, n_orders=n, **conv_kw)
            dR = f"{abs(Ri-prev_R):.1e}" if prev_R is not None else "—"
            dT = f"{abs(Ti-prev_T):.1e}" if prev_T is not None else "—"
            mark = " ◄" if n == N_ORDERS else ""
            print(f"  {n:5d} {2*n+1:6d} {Ri:10.6f} {Ti:10.6f} {Ai:10.6f}"
                  f" {dR:>10} {dT:>10} {ei:10.1e}{mark}")
            prev_R, prev_T = Ri, Ti

        if N_ORDERS >= 45:
            print(f"\n  ✅ n_orders={N_ORDERS} sufficient / 截断阶数足够")
        elif N_ORDERS >= 21:
            print(f"\n  ⚠️  n_orders={N_ORDERS} marginal — consider 45 / 建议提高到45")
        else:
            print(f"\n  ❌ n_orders={N_ORDERS} too low — results may be wrong / 太低，结果可能不准")

    # ── Step 3: Compute / 计算 ──

    # 3A: Specific wavelengths / 指定波长
    results_specific = []
    if WL_MODE in ("A", "C"):
        print_header("Step 3A: Specific Wavelengths / 指定波长计算")
        print(f"\n  Testing {len(WL_SPECIFIC)} wavelengths...")
        print(f"  测试 {len(WL_SPECIFIC)} 个波长...\n")

        for wl in WL_SPECIFIC:
            eps_h = get_si_eps(wl) if EPS_HIGH is None else EPS_HIGH
            kw = dict(common)
            kw['eps_high'] = eps_h
            R, T, A, err = rcwa_te(wl_nm=float(wl), **kw)
            results_specific.append((float(wl), R, T, A, err))

        i_peak = print_table(WL_SPECIFIC, results_specific)
        print(f"\n  Peak A: {results_specific[i_peak][3]*100:.2f}%"
              f" @ {results_specific[i_peak][0]:.1f} nm")

    # 3B: Range scan / 范围扫描
    results_scan = []
    if WL_MODE in ("B", "C"):
        print_header(f"Step 3B: Range Scan {SCAN_START:.0f}–{SCAN_END:.0f} nm / 范围扫描")

        wl_array = np.arange(SCAN_START, SCAN_END + SCAN_STEP / 2, SCAN_STEP)
        print(f"\n  {len(wl_array)} points, step = {SCAN_STEP} nm")
        print(f"  {len(wl_array)} 个点，步长 = {SCAN_STEP} nm\n")

        for wl in wl_array:
            eps_h = get_si_eps(wl) if EPS_HIGH is None else EPS_HIGH
            kw = dict(common)
            kw['eps_high'] = eps_h
            R, T, A, err = rcwa_te(wl_nm=float(wl), **kw)
            results_scan.append((float(wl), R, T, A, err))

        i_peak = print_table(wl_array, results_scan)
        print(f"\n  Peak A: {results_scan[i_peak][3]*100:.2f}%"
              f" @ {results_scan[i_peak][0]:.1f} nm")

    # ── Step 4: Energy check / 能量检查 ──
    all_results = results_specific + results_scan
    if all_results:
        print_header("Step 4: Energy Conservation Check / 能量守恒检查")

        violations = [(wl, R, T) for wl, R, T, A, c in all_results
                      if R + T > 1.0 + 1e-6]
        max_cons = max(c for _, _, _, _, c in all_results)
        max_cons_wl = [wl for wl, _, _, _, c in all_results
                       if c == max_cons][0]

        print(f"\n  Total points / 总计算点:          {len(all_results)}")
        print(f"  Max conservation error / 最大守恒误差: {max_cons:.2e}"
              f" @ {max_cons_wl:.1f} nm")

        if violations:
            print(f"\n  ❌ {len(violations)} points with R+T > 1 (non-physical):")
            print(f"  ❌ {len(violations)} 个点 R+T > 1（非物理）：")
            for wl, R, T in violations:
                print(f"     λ={wl:.1f} nm  R+T={R+T:.6f}")
        elif max_cons < 1e-10:
            print(f"\n  ✅ EXCELLENT: all errors < 1e-10 (machine precision)")
            print(f"  ✅ 优秀：所有误差 < 1e-10（机器精度）")
        elif max_cons < 1e-6:
            print(f"\n  ✅ GOOD: all errors < 1e-6")
            print(f"  ✅ 良好：所有误差 < 1e-6")
        else:
            print(f"\n  ⚠️  MARGINAL: some errors > 1e-6 — increase n_orders")
            print(f"  ⚠️  偏高：部分误差 > 1e-6 — 建议增大截断阶数")

    # ── Step 5: Save CSV / 保存CSV ──
    if SAVE_CSV and all_results:
        data = np.array(all_results)
        np.savetxt(
            CSV_FILENAME, data, delimiter=",",
            header="wavelength_nm,R,T,A,energy_err",
            comments=""
        )
        print(f"\n  📁 Saved / 已保存: {CSV_FILENAME}")

    # ── Summary / 总结 ──
    elapsed = time.time() - t_start

    # Collect peak info from all results
    if all_results:
        A_all = [r[3] for r in all_results]
        i_best = int(np.argmax(A_all))
        best = all_results[i_best]

    print_header("Summary / 总结")
    print(f"""
  Design / 设计:            {DESIGN_NAME}
  Structure / 结构:         Λ={PERIOD:.1f} nm, H={HEIGHT:.1f} nm, FF={FF:.4f}
  Peak A / A峰值:           {best[3]*100:.2f}% @ {best[0]:.1f} nm
  Peak R / R值:             {best[1]*100:.2f}%
  Peak T / T值:             {best[2]*100:.2f}%
  Energy check / 能量检查:   {"✅ PASS" if not violations else "❌ FAIL"}
  Max error / 最大误差:      {max_cons:.2e}
  Convergence / 收敛阶数:    n_orders={N_ORDERS} ({2*N_ORDERS+1} modes)
  Total time / 总耗时:       {elapsed:.2f} s
  CSV: {"saved → " + CSV_FILENAME if SAVE_CSV else "not saved"}

  ────────────────────────────────────────────────────────────────
  CN-627-Cyborg-Lab | Human-AI Collaborative Validation Tool
  CN-627-Cyborg-Lab | 人机协作验证工具
  https://github.com/CN-627-Cyborg-Lab/cn627-rcwa
  ────────────────────────────────────────────────────────────────
""")


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  RUN / 运行                                                               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
if __name__ == "__main__":
    run()