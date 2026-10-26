# 新会话从这里开始

## 先确认什么

1. 运行 `git status --short --branch`，确认当前分支和未提交修改；
2. 运行 `git branch -a` 与 `git tag -n99`，确认功能分支和里程碑；
3. 阅读 `RESEARCH_STATE.md`，不要把 E07-D 当作已完成；
4. 阅读 `docs/DATA_DICTIONARY.md`，确认真实数据列定义和关键文件哈希；
5. 运行 `tools/verify_provenance.ps1`，核对本地原始目录仍与本仓库记录一致。

## 文档导航

- `docs/EXPERIMENT_INDEX.md`：E00～E07C 的入口、证据、状态和缺失项；
- `docs/TIMELINE.md`：按 2026-10-25/26 组织的实验推进；
- `docs/FINDINGS_AND_LIMITATIONS.md`：已验证发现、失败与局限；
- `docs/DATA_DICTIONARY.md`：列定义、参数和阈值；
- `docs/PROVENANCE.md`：来源目录、真实数据路径与 SHA-256；
- `docs/REPRODUCIBILITY.md`：重跑顺序和公平比较规则；
- `docs/ATTITUDE_STATUS.md`：姿态解算启用条件；
- `docs/BRANCHES_AND_TAGS.md`：分支职责和里程碑语义；
- `plans/E07D_TRUE_FVIM_JOINT_OPTIMIZATION.md`：未完成的下一步。

## 当前研究底线

- 真实 TXT 的第 1/2/3 列是转台外/中/内框角位置（deg），不是角速度；第 4 列为 a 通道，第 7 列为 s 通道；第 6 列仅确认是一个加速度测量通道。
- 没有可靠坐标轴、安装矩阵、零位和完整姿态真值前，不报告俯仰/偏航 RMSE。
- 单频 SNR、非目标谐波残差和相对原始测量 RMSE 都是代理指标，不等价于真实角速度或姿态精度。
- 编码器只作为真值/审核参考，不允许进入滤波或参数寻优形成数据泄漏。
- E07-B 已运行 FVIM-style 有限候选演化，但未达到预设训练门槛；E07-C 是固定分解交叉消融，不是新一轮寻优。

## 建议下一动作

先复核 E07-B 候选的同池评价和 E07-C 的 IMF2 阈值跳变，再按 `plans/E07D_TRUE_FVIM_JOINT_OPTIMIZATION.md` 开发连续配对模态权重。开发完成前不得创建 E07-D 完成标签。若继续姿态方向，必须先补齐坐标轴、安装矩阵、同步和加速度轴数证据。

