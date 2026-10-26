# 仓库与数据策略

## 原始材料保护

以下两个目录是本次同步的只读来源，不在其中初始化 Git，也不删除、重命名或覆盖文件：

- `E:\dak_data\Desktop\实验\jyk_Test\DualGyro_FVIM_MVMD_Lab_v0.1`
- `E:\dak_data\Desktop\实验\jyk_Test\DualGyro_FVIM_MVMD_Lab_v0.2_fix`

仓库使用 `archive/` 保存复制快照。复制后的文件允许在后续功能分支中改进，但任何修改必须与原始快照区分。

## 纳入 Git 的内容

- 两版工程中的源代码、README、CSV、TXT、PNG、MAT 和现有结果；
- 合成数据及其 truth MAT；
- E06/E07 已存在的必要中间 MAT，以便复现后续实验；
- 本次整理新增的中文索引、边界说明和只读核验脚本。

## 不纳入 Git 的内容

- 229 个真实 20200824 TXT 原始文件；
- 本机凭据、环境文件和个人配置；
- 未审查的未来批量输出、缓存、临时导出和 MATLAB 自动保存文件。

真实数据由绝对路径和 SHA-256 绑定。需要跨机器复现时，应通过受控数据盘、机构存储或单独授权的数据仓库传递；在明确许可前不使用 Git LFS 上传这些原始 TXT。

## Git LFS

仓库使用 Git LFS 跟踪 `*.mat` 和 `*.png`。克隆后应先确认 `git lfs version`，再执行 `git lfs pull`。CSV、TXT、MATLAB 源码和 Markdown 保持普通 Git 文本对象，便于审阅差异。

