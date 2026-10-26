function cfg = E07A_config()
% E07-A: small-budget physical-retention + noise-proxy fitness feasibility test.
% Compare two candidate alphas ON TRAIN ONLY; do not tune on Z3-2.
root = fileparts(mfilename('fullpath'));
cfg.inputTrain = fullfile(root, 'inputs', 'Z3-1_E06B_record.mat');
cfg.inputTest = fullfile(root, 'inputs', 'Z3-2_E06B_record.mat');
cfg.resultDir = fullfile(root,'results_E07A');
cfg.fs = 1000;
cfg.rollHz = 3;
cfg.mvmd.K = 8;
cfg.mvmd.alphaBaseline = 154.157959862;
cfg.mvmd.alphaCandidates = [80, 300]; % Frozen BEFORE test: exactly two new trials.
cfg.mvmd.maxIterations = 1200;
cfg.mvmd.tolerance = 1e-5;
cfg.mvmd.tau = 0;
cfg.class.csKeep = 0.22;
cfg.class.sampKeep = 1.5;
cfg.class.csMixed = 0.08;
cfg.class.sampMixed = 1.8;
cfg.class.maxSampPoints = 320;
cfg.class.sampTolerance = 0.2;
cfg.wavelet.levels = 4;
cfg.wavelet.thresholdScale = 0.6;
% All predeclared, not modified after seeing held-out results.
cfg.pass.minSidebandRetentionPct = 99;
cfg.pass.maxSidebandRetentionPct = 101;
cfg.pass.maxSidebandPhaseShiftDeg = 0.2;
cfg.pass.minExtraSuppressionDbBothChannels = 1.0;
end
