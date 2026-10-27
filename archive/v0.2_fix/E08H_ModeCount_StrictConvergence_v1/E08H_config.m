function cfg=E08H_config()
% E08-H: developer experiment, NOT automatic FVIM or blind generalization.
root=fileparts(mfilename('fullpath'));
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验'; % 修改成真实TXT根目录
cfg.outputDir=fullfile(root,'results_E08H');
cfg.fs=1000;cfg.rollHz=9;
cfg.trainFiles={'W9-10.TXT','Z9-1.TXT'};
cfg.developmentFiles={'W9-60.TXT','Z9-3.TXT'}; % 已经看过，只能做开发复核
cfg.Kvalues=[3 4 5 6 7 8];
cfg.alphaValues=[154.157959862 320];
cfg.referenceCandidates=[8 154.157959862;8 320];
cfg.maxIter=2400; cfg.tol=1e-5;cfg.tau=0;cfg.eta=1;
cfg.maxAcLoss=.15;cfg.minSamples=2000;
cfg.lastSeconds=3;cfg.skipStartSec=.5;cfg.endFraction=.95;
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
cfg.joint.weights=[.55 .25 .12 .08];
cfg.joint.constraintPenalty=10;cfg.joint.amplitudeScalePct=1;
cfg.joint.phaseScaleDeg=.2;cfg.joint.phaseFloorRawAmpFraction=.01;
cfg.pass.minRetentionPct=99;cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=.2;
end
