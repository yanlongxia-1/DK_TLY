function cfg=E08E_config()
%E08E_CONFIG 只修改dataRoot，然后运行E08E_selftest/E08E_run。
% FVIM-style = 可审计四领先者工程搜索，不冒称原论文完整FVIM源码。
root=fileparts(mfilename('fullpath'));
cfg=E08C_config();
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验'; % 修改为你的TXT上级目录，支持子目录
cfg.outputDir=fullfile(root,'results_E08E_9Hz');
cfg.calibrationCsv=fullfile(root,'calibration_final.csv');
cfg.frequenciesHz=9; % 首轮只做9Hz；第二轮可修改为[3 9 12]
cfg.search.initial=[8 154.157959862;6 185.0015;7 250;8 320]; % 预登记共同种子；不是自动结果表
cfg.search.proposalsPerObjective=4;
cfg.search.maxUniqueCandidates=4+2*cfg.search.proposalsPerObjective;
cfg.search.seed=20261026;
cfg.mvmd.maxIterations=1200; % 维持已有收敛标准，不能在测试集放宽
cfg.mvmd.tolerance=1e-5;
cfg.mvmd.maxAcReconstructionLoss=0.15;
cfg.search.Kmin=5;cfg.search.Kmax=9;
cfg.search.alphaMin=80;cfg.search.alphaMax=450;
cfg.testWspeed=60;         % Wf-60：本次外部角速度评价；训练为Wf-10
cfg.testZmodulation=3;   % Zf-3：若不存在则明确记录缺失，不偷偷换成已看过的Zf-2
cfg.runDevelopmentDiagnostic=false; % true时额外查看Wf-20，不能叫盲测
cfg.trainWindowSeconds=3.2;
cfg.trainTrimStartSec=.5;
cfg.testMaxSeconds=3.0;
cfg.testEvalFraction=[.60 .95];
cfg.minTestSeconds=2.0;
cfg.phase.mode='encoder_initial'; % 台架初始相位仅用于正式角速度评价，搜索目标不接触编码器真值
cfg.huberTune=1.345;cfg.huberMaxIterations=50;cfg.confidenceZ=1.96;
cfg.huberConvergenceTolerance=1e-10;cfg.outlierMadThreshold=6;
cfg.EtrimStartSec=.5;cfg.EtrimEndSec=.2;
cfg.makeFigures=true;
end
