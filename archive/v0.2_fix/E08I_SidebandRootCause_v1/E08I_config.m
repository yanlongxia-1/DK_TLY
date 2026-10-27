function cfg = E08I_config()
% E08-I 预登记：对动态边带失败归因 + 小K严格收敛诊断（开发数据，不是盲测）
root = fileparts(mfilename('fullpath'));
cfg.dataRoot = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验'; % 修改成实际TXT数据总目录
cfg.outputDir = fullfile(root, 'results_E08I');
cfg.fs = 1000; cfg.rollHz = 9; cfg.maxIter = 2400;
cfg.tol = 1e-5; cfg.tau = 0; cfg.eta = 1;
cfg.maxAcLoss = .15; cfg.minSamples = 2000;
cfg.lastSeconds = 3; cfg.skipStartSec = .5; cfg.endFraction = .95;
% 先核查既有 E08-H 保真失败的边界候选，不在开发数据上选择K/alpha
cfg.Zcandidates = [3 154.157959862;4 154.157959862;5 154.157959862;3 320;8 320];
% 原 E08-H 未在 W9-60 上测试小 K。这一步只检查是否严格收敛，绝不据此选参。
cfg.Wcandidates = [3 154.157959862;4 154.157959862;5 154.157959862;6 154.157959862;7 154.157959862;...
                   3 320;4 320;5 320;6 320;7 320];
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
cfg.joint.phaseFloorRawAmpFraction=.01;
cfg.pass.minRetentionPct=99;cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=.2;
end
