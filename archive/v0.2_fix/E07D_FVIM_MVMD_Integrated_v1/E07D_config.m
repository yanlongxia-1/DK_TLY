function cfg=E07D_config()
% E07-D: 联合软模态重构 + 4引导向量FVIM风格K/alpha优化。
% 论文FVIM公开的(7)-(11)更新关系可参考，但作者未公布精确源码/边界操作。
% 本实现明确为可复现工程版本，不能标注为论文作者代码的逐位复现。
root=fileparts(mfilename('fullpath'));
cfg.root=root;
cfg.inputTrain=fullfile(root,'inputs','Z3-1_E06B_record.mat');
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验'; % 填写20200824 TXT根目录才运行跨频率独立检查；留空只运行训练
cfg.outputDir=fullfile(root,'results_E07D');
cfg.fs=1000;cfg.rollHz=3;
cfg.trainFrequenciesHz=[2 3 4];cfg.trainSidebandsHz=[2 4];
% 可先诊断 / 正式优化；已有E07-C证据，默认直接进入有限预算的正式寻优。
cfg.search.Kmin=5;cfg.search.Kmax=9;
cfg.search.alphaMin=80;cfg.search.alphaMax=450;
cfg.search.initialK=[8 6 7 5];
cfg.search.initialAlpha=[154.157959862 115 300 425];
cfg.search.offspringPerMethod=8; % 两目标各8个新候选，总上限4+8+8=20
cfg.search.maxUniqueSolves=20;
cfg.search.randomSeed=20261026;
cfg.mvmd.maxIterations=1200;cfg.mvmd.tolerance=1e-5;cfg.mvmd.tau=0;
cfg.joint.maxAcReconstructionLoss=0.15;
% 共享同一MVMD实现，同一轮候选数据对两个目标函数均评估。
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
% 连续双通道共享权重，非事后按Z3-2调参。阈值中心沿用原方法。
cfg.soft.csWidth=.045;cfg.soft.sampWidth=.085;
cfg.soft.mixExponent=2; % 对接近高熵边界的混叠模态平滑衰减
cfg.joint.weights=[.55 .25 .12 .08];
cfg.joint.constraintPenalty=10; % 动态边带不满足预定保真窗时的显式不可行惩罚
cfg.joint.amplitudeScalePct=1.0;cfg.joint.phaseScaleDeg=.2;
cfg.joint.phaseFloorRawAmpFraction=.01;
cfg.pass.minRetentionPct=99;cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=.2;
% 不根据外部检查文件改搜索参数。本地尚未提供完整TXT时只生成训练结果。
cfg.holdout={'Z6-2.TXT','W6-20.TXT'};
cfg.holdoutRollHz=[6 6];
cfg.holdoutMotionHz=[2 0];
cfg.holdoutTrimStart=.5;cfg.holdoutTrimEnd=.2;
% 公开真实原始数据若已参与此前探索，将在论文注明非盲测。
end
