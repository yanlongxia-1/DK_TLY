function cfg = E06B_config()
%E06B_CONFIG E06-B 一次性联合细节抑噪试验，参数固定、不进行寻优。
root = fileparts(mfilename('fullpath'));
% 改成已有 E06-A results_E06A 文件夹的实际路径（含 Z3-1/Z3-2 子目录）。
cfg.previousResults = 'E:\dak_data\Desktop\实验\jyk_Test\DualGyro_FVIM_MVMD_Lab_v0.2_fix\E06A_Dynamic_Sideband_v1\results_E06A';
cfg.outputRoot = fullfile(root, 'results_E06B');
cfg.cases = {'Z3-1','Z3-2'};
% 只处理已保留的 IMF1，双路共享 Haar 细节系数的收缩系数。
cfg.detailLevels = 6;  % 1000Hz下主要作用于粗略>8Hz的细节，1~5Hz不直接阈值化
cfg.thresholdFactor = 0.75; % 单次先验设置，禁止按 Z3-2 结果重新调参
cfg.minRetentionPct = 99.0; % 两组数据、两路、每个显著动态边带对原始幅值的最低保持率
cfg.maxSidebandPhaseShiftDeg = 0.2; % 候选相对原始边带的绝对相位差上限
cfg.minExtraSuppressionDb = 1.0; % 对原冻结MVMD再提升的非目标谐波残差代理，需两路均达到
end
