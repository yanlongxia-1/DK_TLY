function cfg = E08J_config()
% E08-J v1: physical magnitude/phase fidelity rule audit, no FVIM searching.
% 10% strong/weak boundary is an EXPLORATORY design choice based on E08-I.
% Freeze these numbers BEFORE any independent held-out evaluation.
root = fileparts(mfilename('fullpath'));
cfg.dataRoot = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.outputDir = fullfile(root,'results_E08J');
cfg.fs=1000; cfg.rollHz=9; cfg.maxIter=2400; cfg.tol=1e-5;
cfg.tau=0; cfg.eta=1; cfg.maxAcLoss=0.15;
cfg.minSamples=2000; cfg.lastSeconds=3; cfg.skipStartSec=0.5; cfg.endFraction=0.95;
% Fixed engineering rule; compare hard reconstruction against all-IMF diagnostic.
cfg.class.csKeep=.22; cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08; cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320; cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4; cfg.wavelet.thresholdScale=.6;
% OLD rule maintained for scientific comparison.
cfg.pass.minRetentionPct=99; cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=.2;
% NEW rule: amplitude-reliable components retain old 1% / 0.2 degree limits;
% weak components use a bounded complex-coefficient error in units of the
% strongest raw physical Fourier coefficient IN THE SAME CHANNEL.
cfg.policy.strongFraction=0.10;
cfg.policy.strongAmpTolerancePct=1.0;
cfg.policy.strongPhaseToleranceDeg=0.2;
cfg.policy.strongRelativePhaseToleranceDeg=0.2;
cfg.policy.weakComplexToleranceFraction=0.01;
% Do not use test truth / encoder or the previously viewed W9-60 in training.
cfg.trainFiles={'W9-10.TXT','Z9-1.TXT'};
cfg.Ks=3:8;cfg.alphas=[154.157959862 320];
cfg.enableKnownDevelopment=false; % optional post-freeze audit, NOT blind
cfg.developmentFiles={'W9-60.TXT','Z9-3.TXT'};
end
