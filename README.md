# DK_TLY：双正交陀螺 FVIM–MVMD 研究归档

本仓库用于把本地实验工程、已核实结果、失败记录和下一步计划组织成可追踪的 Git 历史。它不是原论文作者源码，也不把工程代理指标包装成真实姿态精度。

当前同步遵循三条原则：

1. 两个原始实验目录保持只读，仓库仅保存复制后的版本快照；
2. 真实 20200824 TXT 数据不上传，只记录本地路径、文件哈希与复现方法；
3. MAT/PNG 使用 Git LFS，E07-D 未实现前只作为计划存在，不打完成标签。

新会话请先阅读 [START_HERE.md](START_HERE.md)，研究结论与边界见 [RESEARCH_STATE.md](RESEARCH_STATE.md)。实验总表见 [docs/EXPERIMENT_INDEX.md](docs/EXPERIMENT_INDEX.md)，E07-D 真实联合寻优只存在于 [计划](plans/E07D_TRUE_FVIM_JOINT_OPTIMIZATION.md)，尚未实现。

## 目录约定

- `archive/v0.1/`：原始 v0.1 工程快照；
- `archive/v0.2_fix/`：修正版及 E02B～E07C 工程快照；
- `docs/`：实验索引、时间线、数据字典、复现说明和研究边界；
- `tools/`：只读核验脚本；
- `plans/`：尚未完成的研究计划，不代表已有结果。

## 当前结论快照

- v0.1 旧 MVMD 数值更新不收敛，已保留为历史基线；
- v0.2 修复后可在多组参数上收敛，但 E03 物理保护在现有工况未触发；
- E05 只验证特定转台几何下的等效外框输入，固定 KF 在该 RMSE 指标优于冻结 MVMD；
- E06A 显示动态边带约 99%+ 保留，E06B 共享 IMF1 后处理失败；
- E07A/B 均未通过预设训练门槛；
- E07C 证实约 3 dB 代理差主要来自 SampEn=1.8 附近 IMF2 硬分类跳变；
- 下一步是公共候选池、公平双目标评价和连续配对模态权重，而不是把 E07-B 写成成功。

详细数值与限制见 [docs/FINDINGS_AND_LIMITATIONS.md](docs/FINDINGS_AND_LIMITATIONS.md)。

## 数据治理

真实原始数据目录包含 229 个 TXT、约 143 MB，未纳入 Git。关键输入使用 SHA-256 绑定，避免同名文件被误换。仓库内的合成数据、必要中间 MAT 和结果图属于本次研究归档，二进制文件由 Git LFS 管理。

