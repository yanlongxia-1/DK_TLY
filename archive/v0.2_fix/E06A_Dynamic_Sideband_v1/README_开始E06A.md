# E06-A｜真实Z类动态边带保真诊断（冻结双通道MVMD）

这是上一轮计划中提到、现在实际创建的 E06-A MATLAB 工程。

## 实验目的

验证冻结参数的双通道联合MVMD+CS/SampEn模态筛选，在外框存在真实角振动时，是否保留滚转载波 f_roll=3Hz 附近的动态边带 f_roll±f_motion。
本实验只定位动态边带问题，**不修改/优化滤波器**，不复跑229文件标定，不运行FVIM搜索，不对现有KF进行额外对照，也不计算尚未建立坐标轴的俯仰偏航RMSE。

## 运行步骤

1. 解压到独立文件夹，不必覆盖原E05工程。
2. MATLAB打开工程根目录，编辑 `E06_config.m` 的 `cfg.dataRoot`，指定包含`Z3-1.TXT`和`Z3-2.TXT`的20200824原始数据文件夹（可含子目录）。
3. 运行：`E06_dynamic_sideband`
4. 输出将保存在 `results_E06A/Z3-1/` 与 `results_E06A/Z3-2/`。

## 固定参数（源于前期W3-10调参）

- fs=1000 Hz, K=8, alpha=154.157959862, tau=0, 最大1200迭代, tol=1e-5。
- 每组最长截取4秒并联合分解**一次**；零偏仅保持原始数据DC，不参与本次边带判别。
- 复用E05: `gyro_mvmd.m`, `gyro_reconstruct.m`, `gyro_sampen.m`, `gyro_haar.m`, `read_gyro_txt.m`。
- 不收敛时只保存`E06_solver.csv`，不输出失真的最终比较。

## 输出（每个Z工况）

- `E06_solver.csv`: 收敛及负载。
- `E06_mode_details.csv`: 原基准算法模式分类。
- `E06_IMF_modes.csv`: 每个IMF在3Hz和两个候选边带的谐波幅值/相位。
- `E06_IMF_spectrum.png`: 两路真实信号与各IMF频谱。
- `E06_sideband_metrics.csv`: 原始 vs 滤后边带幅值保持率、相位变化。
- `E06_channel_consistency.csv`: 边带及载波处两路相对相位变化。
- `E06_noise_metrics.csv`: 非目标频率残差RMS代理（不是无噪声真值）。
- `E06_spectrum_compare.png`和`E06_waveform_compare.png`: 频谱和时域比较。
- `E06_record.mat`: 保留原始、滤后、分解模态和完整结果，后续E06-B可直接复用而**不用重算MVMD**。
- `E06_interpretation.txt`: 解释范围和风险说明。

## 解读规则与下一步

- `Z3-1`的理论候选边带2和4Hz；`Z3-2`为1和5Hz；实际幅值需与原始及外框编码器运动核对。
- 如果“可靠可辨的真实边带”明显减弱（预设提示阈值＞10%），进入E06-B改进“联合有效边带筛选+软重构”。
- 如果保留良好，直接记录其结论，不强行制造创新；转向双路联合重构的噪声/幅相保真改进。
- **不可把原始信号的单频系数当作无噪声真值**。脚本`heuristicSidebandResolvable`是基于残差白噪声假设的诊断标志，不是严格置信结论。
- 两个工况同样冻结参数，不允许看过Z3-2的结果后改算法仍称Z3-2独立验证；如果需要做创新开发，应另留新的Z工况作为最终盲测。

## 当前程序已检查但尚未MATLAB实机运行

本环境无法运行 MATLAB。若有报错，请发送报错原文、文件名与行号；本次不提供未经实机确认的性能数值。
