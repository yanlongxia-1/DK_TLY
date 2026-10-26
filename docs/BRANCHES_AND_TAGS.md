# 分支与标签

## 分支职责

| 分支 | 提交 | 内容 |
|---|---|---|
| `main` | 持续集成 | 所有已归档成果、状态文档和未完成计划 |
| `feature/baseline-v0.1` | `c0689ba` | v0.1、E02、E02A 历史基线 |
| `feature/v0.2-diagnostics` | `f770ff0` | v0.2 数值修复、E02B/C、E03A-E |
| `feature/calibration-real-data` | `3760818` | E04/CAL01/E05、固定标定与参考姿态接口 |
| `feature/joint-filtering-e06` | `546b8e6` | E06A 动态边带、E06B 联合 IMF1 候选 |
| `feature/fvim-e07abc-diagnostics` | `a14edc9` | E07A/B/C、FVIM-style 与 IMF2 因果消融 |
| `research/attitude-reference-audit` | `3093788` | 姿态解算物理条件与缺口 |
| `plan/e07d-true-fvim-joint-optimization` | `0e422da` | E07-D 设计与验收；不是实现分支 |
| `docs/research-synthesis` | 最终文档提交 | 索引、时间线、复现、来源与状态总结 |

## annotated 标签

| 标签 | 对应完成节点 | 明确不代表什么 |
|---|---|---|
| `milestone-v0.1-baseline` | v0.1 历史基线和 E02A 失败诊断已归档 | 不代表旧 E02 收敛或可作正式基线 |
| `milestone-v0.2-diagnostics` | v0.2 数值修复和 E02/E03 诊断材料完整 | 不代表真实姿态精度提升 |
| `milestone-e05-calibrated-benchmark` | 固定标定与 W3-20 等效输入对照完整 | 不代表任意 body 角速度或姿态解算完成 |
| `milestone-e06-dynamic-diagnostics` | E06A/B 完整运行与停止决策归档 | 不代表 E06B 改进成功 |
| `milestone-e07c-causal-diagnostic` | E07A-C 与 IMF2 机制诊断完整 | 不代表 E07-D 完成或新方法通过盲测 |

所有标签都是 annotated tag，说明文本为中文。计划和文档分支不打“完成”标签。

