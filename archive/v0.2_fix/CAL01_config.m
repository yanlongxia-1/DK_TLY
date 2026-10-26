function cfg = CAL01_config()
% CAL01: per-condition effective bias and gain CANDIDATES, training-only fit.
% No FVIM, MVMD, KF, optimization toolbox, or other dependencies.
root = fileparts(mfilename('fullpath'));
cfg.fs = 1000;
cfg.trainFraction = [0.10, 0.50];
cfg.checkFraction = [0.60, 0.95];
cfg.outputDir = fullfile(root, 'CAL01_results');
% Set real paths below, or place TXT files in the data directory.
cfg.cases = struct('name', {}, 'file', {}, 'fRollHz', {}, 'nominalOuterDegS', {});
cfg.cases(1).name='W3-10';
cfg.cases(1).file=fullfile(root,'data','W3-10.txt');
cfg.cases(1).fRollHz=3.0;
cfg.cases(1).nominalOuterDegS=10;
cfg.cases(2).name='W3-20';
cfg.cases(2).file=fullfile(root,'data','W3-20.txt');
cfg.cases(2).fRollHz=3.0;
cfg.cases(2).nominalOuterDegS=20;
% This flag documents an independently verified geometry hypothesis:
% |Omega_transverse| == |outer gimbal angular rate| when middle frame is fixed.
% DO NOT set true solely because K values appear similar.
cfg.geometryVerified = false;
end
