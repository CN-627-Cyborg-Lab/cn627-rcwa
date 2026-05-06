# CN-627 RCWA Validator
# CN-627 RCWA 验证器

**Developed by CN-627-Cyborg-Lab with AI-assisted code generation (Claude Opus 4.6).**

**由 CN-627-Cyborg-Lab 开发，使用 AI 辅助代码生成（Claude Opus 4.6）。**

Lightweight 1D-RCWA for binary gratings — built for **validation**, not optimization.

轻量级一维RCWA二元光栅计算工具 — 专为**验证**而非优化设计。

---

## Features / 特性

- **TE polarization, normal incidence** — clean, auditable implementation
  **TE偏振，正入射** — 清晰可审计的实现
- **TMM cross-check** — flat slab (FF=1) matches TMM to machine precision (~1e-15)
  **TMM交叉验证** — 平板结构（FF=1）与TMM一致至机器精度（~1e-15）
- **Energy conservation** — |R+T-1| < 1e-10 for lossless structures
  **能量守恒** — 无损耗结构 |R+T-1| < 1e-10
- **Pure NumPy** — no hidden dependencies
  **纯NumPy** — 无隐藏依赖
- **Si dispersion** — Green 2008 data built-in (250–1200 nm)
  **硅色散** — 内置Green 2008数据（250-1200 nm）
---

## Development Model / 开发模式

This project was developed through human-AI collaboration.

本项目通过人机协作模式开发。

| Role / 角色 | Human (T.z) / 人类 | AI (Claude Opus 4.6) |
|------|-------------|----------------------|
| Physics / 物理 | Method selection, validation / 方法选择、验证 | No autonomous physics claims / 无自主物理主张 |
| Code / 代码 | No manual coding / 无手动编码 | Full implementation / 完整实现 |
| Debugging / 调试 | Result verification, direction / 结果验证、方向指导 | Code generation, error fixing / 代码生成、错误修复 |
| Verification / 验证 | TMM cross-check, energy conservation / TMM对比、能量守恒 | Test script generation / 测试脚本生成 |

All AI-generated code was validated against analytical solutions (TMM) and published literature (Moharam 1995, Li 1996).

所有 AI 生成代码均经解析解（TMM）和已发表文献（Moharam 1995, Li 1996）验证。

The author's contribution lies in scientific supervision, validation, and direction, not manual coding.

作者的贡献在于科学监督、验证和方向指导，而非手动编码。

---

## Install / 安装
```bash
pip install cn627-rcwa
```

Or from source / 或从源码安装:
```bash
git clone https://github.com/CN-627-Cyborg-Lab/cn627-rcwa.git
cd cn627-rcwa
pip install -e ".[test]"
```

---

## 5-Minute Quick Start / 5分钟快速入门

### Option A: Zero-code validation / 零代码验证

The fastest way — edit parameters, run, read report:

最快方式 — 改参数，运行，看报告：
```bash
# 1. Open and edit your parameters / 打开并修改参数
#    examples/00_validate_your_design.py  →  Section 1 & 2

# 2. Run / 运行
python examples/00_validate_your_design.py

# 3. Read the report in terminal / 在终端阅读报告
```

What you need to edit / 你需要修改的：
```python
# Section 1: Your structure / 你的结构
PERIOD = 500.0                          # nm, grating period / 光栅周期
HEIGHT = 600.0                          # nm, grating height / 光栅高度
FF = 0.50                               # fill fraction / 占空比

# Materials: n,k → eps = complex(n,k)**2 / 材料转换
EPS_HIGH = complex(3.48, 0.00548) ** 2   # high-index / 高折射率
EPS_LOW  = 1.0                           # low-index (Air) / 低折射率（空气）
EPS_INC  = 1.0                           # incidence (Air) / 入射（空气）
EPS_SUB  = 1.47 ** 2                     # substrate (SiO2) / 基底

# Section 2: Which wavelengths / 测哪些波长
WL_MODE = "C"                            # "A"=specific, "B"=scan, "C"=both
WL_SPECIFIC = [850, 905, 950, 1000]      # specific points / 指定波长
SCAN_START = 900.0                       # scan range / 扫描范围
SCAN_END = 910.0
SCAN_STEP = 0.2
```

The script will automatically / 脚本会自动：
- ✅ Verify convergence / 验证收敛性
- ✅ Check energy conservation / 检查能量守恒
- ✅ Flag non-physical results / 标记非物理结果
- ✅ Save CSV / 保存CSV

### Option B: Python API / Python接口
```python
from cn627 import rcwa_te
import numpy as np

# Define material / 定义材料
n, k = 3.48, 0.00548
eps_si = complex(n, k) ** 2    # n,k → permittivity / 折射率→介电常数

# Single wavelength / 单波长
R, T, A, err = rcwa_te(
    wl_nm=905.0,
    period_nm=500.0,
    height_nm=600.0,
    ff=0.50,
    eps_high=eps_si,
    eps_sub=1.47**2,
    n_orders=45
)
print(f"R={R:.4f}  T={T:.4f}  A={A:.4f}  err={err:.1e}")
```
```python
# Spectral scan / 光谱扫描
from cn627 import rcwa_spectrum

wl = np.arange(900, 911, 0.2)
result = rcwa_spectrum(wl, period_nm=500.0, height_nm=600.0, ff=0.50,
                       eps_high=eps_si, eps_sub=1.47**2, n_orders=45)

i = np.argmax(result['A'])
print(f"Peak A: {result['A'][i]*100:.1f}% @ {result['wl'][i]:.1f} nm")
```

### Material conversion cheat sheet / 材料转换速查表

| Material / 材料 | n | k | eps = complex(n,k)**2 |
|------|---|---|---|
| Air / 空气 | 1.00 | 0 | `1.0` |
| SiO2 | 1.47 | 0 | `1.47**2` |
| Si (905nm) | 3.48 | 0.00548 | `complex(3.48,0.00548)**2` |
| Ge | 4.00 | 0 | `16.0` |
| Se | 2.64 | 0 | `6.9696` |
| Your material / 你的材料 | n | k | `complex(n,k)**2` |

---
### Structure diagram / 结构示意图

```text
┌─────────────────────────────┐
│  eps_inc (Air = 1.0)        │  ← light in / 光入射
├─────────────┬───────────────┤
│ HIGH │ LOW  │ HIGH │ LOW    │  ← grating / 光栅
│      │      │      │        │     Λ, H, FF
├──────┴──────┴──────┴────────┤
│  eps_sub (SiO2 = 1.47²)     │  ← light out / 光透射
└─────────────────────────────┘
```
---

## Validation / 验证
```bash
pytest tests/ -v
```
---

| Test / 测试 | What it checks / 检验内容 |
|------|----------------|
| `test_flat_slab` | RCWA(FF=1) vs TMM, multiple materials / 多种材料平板对比 |
| `test_energy` | R+T ≤ 1 (absorbing); R+T = 1 (lossless) / 能量守恒 |
| `test_convergence` | Results stabilize with increasing orders / 截断阶数收敛 |

---

## Examples / 示例

| Script / 脚本 | Description / 描述 |
|--------|-------------|
| `00_validate_your_design.py` | One-stop validation tool / 一站式验证工具 |
| `01_si_grating.py` | Si/Air grating, visible range / 硅/空气光栅，可见光范围 |
| `02_device1_validation.py` | Niraula 2014 Device 1 (Ge/Se @ 10.6 µm) / 论文验证 |

---

## Algorithm / 算法

S-matrix formulation: Moharam et al., JOSA A 12(5), 1068 (1995).
S矩阵公式：Moharam 等，JOSA A 12(5), 1068 (1995)。

Redheffer star product: Li, JOSA A 13(5), 1024 (1996).
Redheffer星积：Li，JOSA A 13(5), 1024 (1996)。

---

## Limitations / 限制

- TE only (TM planned for v0.2) / 仅TE偏振（TM计划于v0.2）
- Normal incidence only / 仅正入射
- Single grating layer / 单层光栅
- 1D periodicity / 一维周期

---

## Citation / 引用
```bibtex
@software{cn627_rcwa,
  author  = {T.z},
  title   = {CN-627 RCWA Validator},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/CN-627-Cyborg-Lab/cn627-rcwa}
}
```

---

## License / 许可证

MIT