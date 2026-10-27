function cfg=E08K_config()
% E08-K v1: 9 Hz FVIM-STYLE paired search and fixed-calibration rate benchmark.
% Explicitly NOT the original FVIM author's exact implementation.
% Search never reads encoder truth or test data. Freeze BEFORE validation.
root=fileparts(mfilename('fullpath'));
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.outputDir=fullfile(root,'results_E08K_9Hz');
cfg.calibrationCsv=fullfile(root,'calibration_final.csv');
cfg.fs=1000;cfg.frequenciesHz=9;
cfg.search.initial=[8 154.157959862; 5 154.157959862; 7 450; 8 320];
cfg.search.proposalsPerObjective=4;
cfg.search.maxUniqueCandidates=4+2*cfg.search.proposalsPerObjective;
cfg.search.seed=20261027;
cfg.search.Kmin=3;cfg.search.Kmax=8;
cfg.search.alphaMin=80;cfg.search.alphaMax=450;
cfg.mvmd.maxIterations=2400;cfg.mvmd.tolerance=1e-5;
cfg.mvmd.tau=0;cfg.mvmd.maxAcReconstructionLoss=.15;
% Same E08-J hard reconstruction: no rule changes in this experiment.
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
cfg.pass.minRetentionPct=99;cfg.pass.maxRetentionPct=101;
cfg.pass.maxSidebandPhaseShiftDeg=.2;
% Fixed exploratory physical-fidelity rule from E08-J; no tuning on W9-60 / Z9-3.
cfg.policy.strongFraction=.10;
cfg.policy.strongAmpTolerancePct=1.0;
cfg.policy.strongPhaseToleranceDeg=.2;
cfg.policy.strongRelativePhaseToleranceDeg=.2;
cfg.policy.weakComplexToleranceFraction=.01;
% Both train and tests use the same predetermined E08-J 3000-sample tail rule.
cfg.lastSeconds=3;cfg.skipStartSec=.5;cfg.endFraction=.95;
cfg.minTestSeconds=2.0;
cfg.testWspeed=60;cfg.testZmodulation=3;
cfg.runDevelopmentDiagnostic=false;
cfg.phase.mode='encoder_initial';
cfg.EtrimStartSec=.5;cfg.EtrimEndSec=.2;
cfg.huberTune=1.345;cfg.huberMaxIterations=50;
cfg.huberConvergenceTolerance=1e-10;
cfg.outlierMadThreshold=6;cfg.confidenceZ=1.96;
cfg.makeFigures=true;
end
