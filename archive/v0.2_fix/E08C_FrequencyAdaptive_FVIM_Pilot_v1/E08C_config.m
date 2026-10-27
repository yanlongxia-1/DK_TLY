function cfg = E08C_config()
% E08-C：跨滚转频率参数适应性可行性实验（不是全频率最终验证）。
% 仅需修改 dataRoot。不得在验证工况上反向调节参数或权重。
root=fileparts(mfilename('fullpath'));
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.outputDir=fullfile(root,'results_E08C');
cfg.fs=1000;
cfg.frequenciesHz=[3 9]; % 先研究差异明显的3和9Hz；验证后再扩展6/12/21Hz
cfg.trainWsuffix='-10.TXT';
cfg.trainZsuffix='-1.TXT';
cfg.checkWsuffix='-20.TXT';
cfg.checkZsuffix='-2.TXT';
cfg.trimStartSec=.5;
cfg.trimEndSec=.2;
cfg.maxWindowSeconds=3.2; % 防止大规模FVIM占用内存
cfg.mvmd.maxIterations=1200;
cfg.mvmd.tolerance=1e-5;
cfg.mvmd.tau=0;
cfg.mvmd.maxAcReconstructionLoss=.15;
cfg.search.initial=[8 154.157959862; 6 185.0015; 7 250; 8 320];
cfg.search.Kmin=5;cfg.search.Kmax=9;
cfg.search.alphaMin=80;cfg.search.alphaMax=450;
cfg.search.proposals=4;  % 每频率4个初始+4个更新，共8个候选
cfg.search.seed=20261026;
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
cfg.joint.weights=[.55 .25 .12 .08];
cfg.joint.constraintPenalty=10;
cfg.joint.amplitudeScalePct=1;
cfg.joint.phaseScaleDeg=.2;
cfg.joint.phaseFloorRawAmpFraction=.01;
cfg.pass.minRetentionPct=99;cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=.2;
end
