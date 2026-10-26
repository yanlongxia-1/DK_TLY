function cfg = E06_config()
% E06-A: Z类动态边带检验配置。只修改 dataRoot 即可。
root = fileparts(mfilename('fullpath'));
cfg.dataRoot = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';  % <-- 修改为包含 Z3-1.TXT、Z3-2.TXT 的目录
cfg.outputRoot = fullfile(root, 'results_E06A');
cfg.testFiles = {'Z3-1.TXT', 'Z3-2.TXT'};
cfg.fs = 1000;                  % 名义采样率；历史档案的1000Hz估计
cfg.rollHz = 3;
cfg.skipSeconds = 0.8;          % 去掉起始过渡段
cfg.maxDurationSeconds = 4.0;  % 最多4秒/文件，降低MVMD运行负担
cfg.minDurationSeconds = 2.5;
% 完全冻结 E05/W3-10 的基准配置，本次不执行 FVIM。
cfg.mvmd.K = 8;
cfg.mvmd.alpha = 154.1579598620;
cfg.mvmd.maxIterations = 1200;
cfg.mvmd.tolerance = 1e-5;
cfg.mvmd.tau = 0;
cfg.class.csKeep = 0.22;
cfg.class.sampKeep = 1.5;
cfg.class.csMixed = 0.08;
cfg.class.sampMixed = 1.8;
cfg.class.maxSampPoints = 320;
cfg.class.sampTolerance = 0.2;
cfg.wavelet.levels = 4;
cfg.wavelet.thresholdScale = 0.6;
% 6dB阈值是分析性标志，不是统计学真实信号证明；不用于调整滤波器。
cfg.sidebandReliabilityDb = 6;
cfg.attenuationAlertPct = 10;
cfg.plotMaxHz = 15;
end
