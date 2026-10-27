function cfg=E08A_config()
%E08A_CONFIG E08-A: 滚转频率/标定参数固定，真实W类外框等效角速度评价。
% 唯一必改项: cfg.dataRoot 指向包含 E6.TXT/E9.TXT/W6-20.TXT/W9-20.TXT 的上级目录。
root=fileparts(mfilename('fullpath'));
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.outputDir=fullfile(root,'results_E08A');
cfg.calibrationCsv=fullfile(root,'calibration_final.csv');
cfg.lockedParamsCsv=fullfile(root,'E07D_locked_params.csv');
cfg.cases = {'W9-20.TXT'};
cfg.mvmd.maxIterations = 2400;
cfg.outputDir = fullfile( ...
    fileparts(mfilename('fullpath')), ...
    'results_E08B_W9');
cfg.fs=1000;
cfg.evalFraction=[0.60 0.95]; % 不在评价区间调参，所有方法同段输入
cfg.EtrimStartSec=0.5;cfg.EtrimEndSec=0.2;
% 与原标定阶段一致的Huber稳健谐波拟合设置。
cfg.huberTune=1.345;cfg.huberMaxIterations=50;
cfg.huberConvergenceTolerance=1e-10;cfg.outlierMadThreshold=6;
cfg.confidenceZ=1.96;
cfg.minimumEvalSec=2.0;
% 滚转角: 仅初始化相位+固定频率积分；不跟踪真实滚转角。
cfg.phase.mode='encoder_initial';  % 用编码器第一样本作为台架初相位，仅离线公平比较
cfg.phase.manualPhi0Deg=NaN;     % 改为'manual'时必须设置真实已知的初相位
cfg.phase.manualRollSign=NaN;    % +1或-1，必须与安装角约定一致
% 仅使用E类纯滚转固定零输入偏置；W工况直流截距差作为限制记录。
% 基线和新算法全部冻结，不允许在W6/W9重寻优。
cfg.methods.entropyPoolDiagnostic=true;
cfg.mvmd.maxIterations=1200;cfg.mvmd.tolerance=1e-5;cfg.mvmd.tau=0;
cfg.mvmd.maxAcReconstructionLoss=0.15;
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
cfg.soft.csWidth=.045;cfg.soft.sampWidth=.085;
cfg.soft.mixExponent=2;
cfg.figureVisible='off';
end
