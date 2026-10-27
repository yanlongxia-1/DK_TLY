function cfg = E08F_config()
% E08F: numerical diagnostics only, NOT testing new K/alpha, NOT a filter comparison.
root = fileparts(mfilename('fullpath'));
cfg.dataRoot = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.outputDir = fullfile(root,'results_E08F');
cfg.fs=1000;
cfg.rollHz=9;
cfg.filenames={'W9-60.TXT','Z9-3.TXT'};
cfg.methods = {'fixed_baseline','entropy_auto','physics_auto'};
cfg.candidates=[8 154.157959862; 7 450; 8 320];  % LOCKED from E08E
cfg.maxIterDiagnostic=2400;     % Diagnostic only: never relax tolerances on held-out data.
cfg.tol=1e-5;
cfg.tau=0;
cfg.maxAcLoss=.15;
cfg.testFraction=[.60 .95];
cfg.testMaxSec=3.0;
cfg.minTestSec=2.0;
cfg.skipStartSec=.5;
cfg.logIterations=[100 200 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400];
end
