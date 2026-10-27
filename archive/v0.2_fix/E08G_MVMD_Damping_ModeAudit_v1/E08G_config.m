function cfg=E08G_config()
% E08-G: FIXED-param solver experiment, not FVIM optimization / held-out evaluation.
root=fileparts(mfilename('fullpath'));
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.outputDir=fullfile(root,'results_E08G');
cfg.fs=1000;cfg.rollHz=9;
cfg.filenames={'W9-60.TXT','Z9-3.TXT'};
cfg.names={'fixed_baseline','entropy_auto','physics_auto'};
cfg.candidates=[8 154.157959862;7 450;8 320];
cfg.relaxValues=[1.0 0.5]; % 1=original; .5=damped-center update (one prespecified intervention)
cfg.maxIter=2400;cfg.tol=1e-5;cfg.tau=0;cfg.maxAcLoss=.15;
cfg.skipStartSec=.5;cfg.endFraction=.95;cfg.lastSeconds=3;
cfg.checkpoints=[600 1200 1800 2400];
% E08E baseline hard reconstruction defaults, unchanged
cfg.class.csKeep=.22;cfg.class.sampKeep=1.5;
cfg.class.csMixed=.08;cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;cfg.class.sampTolerance=.2;
cfg.wavelet.levels=4;cfg.wavelet.thresholdScale=.6;
cfg.physics.bandHz=5;cfg.physics.protectHz=1;cfg.physics.minCarrierShare=.05;cfg.physics.minBandShare=.05;
end
