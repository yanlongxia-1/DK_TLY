function cfg = E07B_config()
% E07-B: 双通道保真 FVIM-style 目标函数；本实验不是原论文作者源码复刻。
root=fileparts(mfilename('fullpath'));
cfg.inputTrain=fullfile(root,'inputs','Z3-1_E06B_record.mat');
cfg.inputCheck=fullfile(root,'inputs','Z3-2_E06B_record.mat');
cfg.outputDir=fullfile(root,'results_E07B');
cfg.fs=1000; cfg.rollHz=3;
% 仅 Z3-1 作为训练，Z3-2 在冻结参数后才可检查，禁止反向调整。
cfg.trainComponentsHz=[2 3 4]; cfg.trainSidebandHz=[2 4];
cfg.checkComponentsHz=[1 3 5]; cfg.checkSidebandHz=[1 5];
% 与已验证算法同一求解器与重构流程。保守限制实际新增运行量。
cfg.mvmd.alphaBaseline=154.157959862;
cfg.mvmd.maxIterations=1200; cfg.mvmd.tolerance=1e-5; cfg.mvmd.tau=0;
cfg.search.Kmin=4; cfg.search.Kmax=8;
cfg.search.alphaMin=110; cfg.search.alphaMax=350;
cfg.search.initialK=[8 6 6 4];
cfg.search.initialAlpha=[154.157959862 110 300 220];
cfg.search.offSpringPerMethod=4;
cfg.search.maxUniqueTrainSolves=12; % 含四个共同初值，每个适应度分支再提出四个
cfg.search.randomSeed=20261026;
% CS/SampEn 工程版模态分类阈值。两组方法完全一致。
cfg.class.csKeep=0.22; cfg.class.sampKeep=1.5;
cfg.class.csMixed=0.08; cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320; cfg.class.sampTolerance=0.2;
cfg.wavelet.levels=4; cfg.wavelet.thresholdScale=0.6;
% 联合物理保真项；完全固定，禁止根据检查集调节。
% 这里的 noiseResidual 是谐波投影的非目标残差代理，绝非真值SNR。
cfg.joint.weights=[0.55 0.25 0.12 0.08];  % 噪声、幅值、相位、相对相位
cfg.joint.amplitudeScalePct=1.0;    % 1%幅值变化作归一化单位
cfg.joint.phaseScaleDeg=0.2;       % 0.2度相位变化作归一化单位
cfg.joint.phaseFloorRawAmpFraction=0.01; % 极小幅值分量不施相位约束
cfg.joint.maxAcReconstructionLoss=0.15;
% 固定验收规则（不是适应度）：阈值来自前轮E07-A方案。
cfg.pass.minRetentionPct=99; cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=0.2;
cfg.pass.minExtraSuppressionDb=1.0; % 两路分别都需≥1dB
cfg.pass.minAgainstEntropyDb=0.3; % 共同合格时比旧适应度至少有可见净改善
cfg.checkOnTrainPassOnly=true;
end
