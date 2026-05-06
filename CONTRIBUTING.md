# Contributing to CN-627 / 贡献指南

## AI-Assisted Development / AI 辅助开发

This project demonstrates that rigorous scientific software can be built through AI-human collaboration, even without formal programming training.

本项目证明严谨的科学软件可以通过人机协作构建，即使没有专业编程训练。

**Model / 模式:**
- AI generates code / AI 生成代码
- Human verifies physics / 人类验证物理
- Transparent disclosure / 透明披露

---

## Bug Reports / 报告Bug

Open an issue with:
请提交issue，包含以下信息：

- Python version, NumPy version / Python版本、NumPy版本
- Minimal reproduction code / 最小复现代码
- Expected vs actual output / 预期输出 vs 实际输出

---

## Code Contributions / 代码贡献

Contributors are welcome to use AI tools, but must:

欢迎贡献者使用 AI 工具，但必须：

1. **Verify against literature / 对照文献验证**
   - All AI-generated code must be validated against analytical solutions (TMM) or published literature.
   - 所有 AI 生成代码必须对照解析解（TMM）或已发表文献验证。

2. **Disclose AI use / 披露 AI 使用**
   - Mention AI tools in your PR description (e.g., "Code generated with Claude, verified against Moharam 1995").
   - 在 PR 描述中提及 AI 工具（例如"使用 Claude 生成，对照 Moharam 1995 验证"）。

3. **Pass all tests / 通过所有测试**
   - Run `pytest tests/ -v` before submitting.
   - Ensure energy conservation: `|R+T-1| < 1e-10` for lossless structures.
   - 提交前运行 `pytest tests/ -v`。
   - 确保能量守恒：无损耗结构 `|R+T-1| < 1e-10`。

4. **Understand your code / 理解你的代码**
   - You must be able to explain every line you submit.
   - AI-generated code you cannot explain will not be merged.
   - 你必须能够解释你提交的每一行代码。
   - 你无法解释的 AI 生成代码不会被合并。

---

## Scope / 范围

Welcome / 欢迎:
- ✅ TM polarization / TM偏振
- ✅ Oblique incidence / 斜入射
- ✅ Multi-layer support / 多层支持
- ✅ Bug fixes and tests / Bug修复和测试
- ✅ Documentation improvements / 文档改进

Out of scope / 不在范围内:
- ❌ Optimization / 优化算法
- ❌ GUI / 图形界面
- ❌ Non-NumPy backends / 非NumPy后端