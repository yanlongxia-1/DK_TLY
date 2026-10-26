# 复现实验说明

## 1. 获取仓库

```powershell
git clone https://github.com/yanlongxia-1/DK_TLY.git
Set-Location DK_TLY
git lfs pull
git status --short --branch
```

如果只核对已归档节点，可检出相应 annotated 标签；标签与分支对应关系见 `docs/BRANCHES_AND_TAGS.md`。

## 2. 核对数据

在原 Windows 机器运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\verify_provenance.ps1
```

脚本只读取两个实验来源目录和 W3-10、W3-20、Z3-1、Z3-2 四个关键真实 TXT，不修改原始数据。

跨机器时，把 20200824 数据通过受控渠道放到本地，再修改各配置中的 `cfg.file`、`cfg.dataRoot` 或 `cfg.previousResults`。不要把个人绝对路径提交回仓库；推荐复制配置为被 `.gitignore` 排除的 `local_config_*.m`，或在独立工作分支修改。

## 3. MATLAB 依赖

代码使用 MATLAB 的 `readmatrix`、table、基础绘图和线性代数；MVMD、Haar、SampEn 等核心实现已自包含。仓库没有记录准确 MATLAB 版本和工具箱清单，因此首次复现时必须在日志中补充 `version`、`ver`、操作系统、CPU、内存和运行时间。

## 4. 推荐顺序

1. v0.2 `E00_selftest.m`；若失败，停止。
2. `E02B_solver_audit.m` 的受控数值审计；不要直接扩大为高负载网格。
3. `E02C_imf_threshold_audit.m`，只用训练段。
4. E03B/C 合成真值，用于确认误差计算和不收敛处理。
5. E04/CAL01/E05；确认原始路径、固定标定表和坐标边界。
6. E06A；已有 MAT 可直接支持 E06B，不必重复分解。
7. E06B、E07A、E07B、E07C 按各自 README；不得在失败后反复调参并仍称预注册验证。
8. E07-D 目前只有计划，先建立新目录和测试，不覆盖 E07A-C。

E07-B 的候选并集同池复核可直接运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\review_e07b_same_pool.ps1
```

## 5. 关键实验入口

- v0.1 历史：`archive/v0.1/`
- v0.2 数值与诊断：`archive/v0.2_fix/`
- E05 标定 W 基准：`archive/v0.2_fix/E05_Calibrated_W_Benchmark_v1/`
- E06A/B：`archive/v0.2_fix/E06A_Dynamic_Sideband_v1/`、`E06B_Joint_IMF1_Pilot_v1/`
- E07A/B/C：`archive/v0.2_fix/E07A_JointObjective_Pilot_v1/`、`E07B_PhysicsFVIM_MVMD_v1/`、`E07C_IMF2_CausalAblation_v1/`

## 6. 公平比较规则

- 使用相同原始文件、相同时间窗、相同固定标定和相同真值；
- 编码器只作结果评价，不能参与滤波、模态选择或优化；
- 分开报告：求解器收敛、代理指标、合成真值误差、真实等效输入误差、姿态误差；
- 报告全部不收敛和失败候选，不只挑最好结果；
- 改过阈值、目标或搜索范围就建立新版本，不覆盖旧结果；
- 最终结论必须包含 raw、冻结基线、改进方法和消融。

## 7. 同步时的验证范围

本次同步完成了文件计数、大小、SHA-256、Git LFS 指针、分支/标签和结果表结构核对。MATLAB 中仅运行并通过了两个轻量自检：`E07B_selftest`（不跑 MVMD）和 `E06B_selftest`（Haar 零阈值可逆性）；没有重跑任何 MVMD/FVIM 或完整性能实验。未来重跑必须把新日志与原结果分目录保存，并在实验索引中标明“历史输出”或“本次重跑”。

