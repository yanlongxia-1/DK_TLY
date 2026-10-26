# 实验索引

本表只登记磁盘中实际存在的代码与输出。`完成`表示该实验节点已有脚本和结果文件，不表示方法优于基线。

| 节点 | 数据性质 | 实际入口/目录 | 已有证据 | 状态与结论 |
|---|---|---|---|---|
| E00 | 合成 | `archive/v0.1`、`archive/v0.2_fix` 的 `E00_selftest.m` | selftest TXT/CSV/PNG | 完成自检材料归档；用于公式与求解器检查 |
| E01 | 真实 W3-10 | v0.1 `E01_data_audit.m` | 数据审计 CSV/频谱图 | 完成历史数据审计；不含真实姿态真值 |
| E02 v0.1 | 真实 W3-10 | `E02_paper_baseline.m` | E02 CSV/MAT/PNG | 历史基线；E02A 表明旧 `tau=0.025...0.25` 六组均未收敛，不作正式论文基线 |
| E02A | 真实 W3-10 | v0.1 `E02A_convergence_diagnostic.m` | convergence/modes/诊断图 | 完成失败诊断；保留而不修饰 |
| E02B_SAFE | 真实 W3-10 训练段 | v0.2 `E02B_solver_audit.m` 对应安全输出 | 9 组 K/alpha 数值诊断 | 9 组中 6 组在 250 次上限内收敛；是安全数值诊断，不是完整 FVIM 证明 |
| E02C | 真实 W3-10 训练段 | `E02C_imf_threshold_audit.m` | 阈值表、模态表、说明 | 完成工程阈值敏感性审计；禁止用验证段反调阈值 |
| E03A | 已查看的 W3-10 验证段 | `E03A_fixed_MVMD_physics_check.m` | summary、mode changes、图 | 物理保护触发 0、分类变化 0、输出差 0；该工况没有增益 |
| E03B | 合成 truth | `E03B_synthetic_truth_test.m` | 三场景真实 RMSE | 冻结滤波相对 raw 降低合成真值 RMSE，但 E02/E03 相同，不能证明物理保护有效 |
| E03C | 合成 truth | `E03C_sideband_audit.m` | 频率×种子诊断表 | 多数场景收敛并降低合成 RMSE；8 Hz 某种子不收敛；E02/E03 仍相同 |
| E03D | 合成 truth | `E03D_phase_calibration_audit.m` | 相位/标定结果与说明 | 完成相位、偏置和增益误差归因研究 |
| E03E | 合成/假设观测 | `E03E_phase_reference_and_calibration.m` | phase seed/summary/标定分离 | 完成需求研究；第 6 列加速度能否提供相位仍未证实 |
| E04 | 真实 W3-20 | `E04_real_reference_audit.m` | 列统计、窗口、候选、图 | 未发现可直接认定为 3 Hz body roll-rate 的转台列；不构造姿态真值 |
| CAL01 | 真实 W3-10/W3-20 | `CAL01_run.m` | 当前快照仅含代码；历史配置指向真实目录 | 有效偏置/增益候选，不等于独立实验室绝对标定 |
| E05 固定参数 | 真实 W3-10/W3-20 | `E05_RealData_FixedCalibration` | 代码、说明；部分运行需外部 TXT | 固定标定缺失时不制造 deg/s；姿态参考默认关闭 |
| E05 标定基准 | 真实 W3-20 | `E05_Calibrated_W_Benchmark_v1` | calibration table、CSV/MAT/PNG | 完成等效外框恒速输入对照；KF RMSE 5.2432、MVMD 5.5913、raw 5.9568 deg/s |
| E06A | 真实 Z3-1/Z3-2 | `E06A_Dynamic_Sideband_v1` | 两工况完整 CSV/MAT/PNG | 两次 MVMD 均收敛；冻结算法保留动态边带约 99%+，原“硬保护误删边带”假设不成立 |
| E06B | 真实 Z3-1/Z3-2 已保存中间量 | `E06B_Joint_IMF1_Pilot_v1` | quality/sidebands/decision/MAT/PNG | 侧带门槛和额外抑噪门槛都未通过，按预注册规则停止 |
| E07A | Z3-1 训练；Z3-2 未运行 | `E07A_JointObjective_Pilot_v1` | solver/objective/decision | 固定 alpha 候选未达到两路各 1 dB 改善，停止，不看 Z3-2 反调 |
| E07B | Z3-1 训练 | `E07B_PhysicsFVIM_MVMD_v1` | 12 候选、两目标、轨迹、求解器、MAT/PNG | FVIM-style 有限演化已运行，但新目标未通过训练条件；Z3-2 未作为盲测运行 |
| E07C | Z3-1 训练 | `E07C_IMF2_CausalAblation_v1` | 2 次分解×2 策略交叉消融 | 完成机制诊断：约 3 dB 代理差主要来自 IMF2 硬分类跳变 |
| E07-D | 尚未锁定的新记录 | `plans/E07D_TRUE_FVIM_JOINT_OPTIMIZATION.md` | 只有计划 | **未实现、未运行、无完成标签** |
| E08～E10 | 合成/静态计划入口 | 根目录脚本仅为薄入口，当前无对应完整新输出 | 少量入口代码 | 不标记为已完成研究节点，缺失项待补 |

## 缺失项

- 原论文作者的完整 FVIM/MVMD 可执行源码和唯一适应度细节；
- CAL01 当前运行输出目录（快照内未发现）；
- E05 固定参数版对 W3-10/W3-20 的完整结果目录；
- E07B `E07B_locked_check.csv`（训练未通过，所以按设计未生成）；
- 新鲜盲测 W/Z 记录的锁定清单与结果；
- E07-D 源码、配置、结果与标签；
- 完整姿态真值、三轴加速度确认和坐标安装矩阵。

