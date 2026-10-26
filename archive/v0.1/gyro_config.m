function cfg = gyro_config()
% Dual Orthogonal MEMS Gyro MVMD Laboratory - unified config
% Only change cfg.file for one-file E01-E03/E05-E06 experiments.
root = fileparts(mfilename('fullpath'));
cfg.root = root;
cfg.file = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验\20200824-1-2\W3-10.txt';
cfg.outputRoot = fullfile(root, 'results');
cfg.fs = 1000;                        % Hz; do NOT guess actual sampling frequency
cfg.cols = [4, 7];                   % historical TXT: col4=a, col7=s
cfg.rollHz = 3;                      % reference only; check in E01
cfg.phi0 = 0;                        % assumed for DEMOD experiments
cfg.biasMode = 'none';                % 'none' or 'manual'; not automatically mean removal
cfg.bias = [-0.186041, -0.072559];   % EXAMPLE ONLY for a/s; verify per dataset
cfg.sensitivity = [NaN, NaN];        % [Ka, Ks]: physical units per deg/s; configure if known
cfg.seed = 20261025;
cfg.staticData = false;                % E10: explicitly verify platform is stationary
cfg.minAllanSeconds = 120;             % 10+ minutes is preferred for useful low-frequency terms

% Temporal split: disjoint train 10%-50%, validation 60%-95%.
cfg.trainFraction = [0.10, 0.50];
cfg.valFraction = [0.60, 0.95];

% MVMD / FVIM parameters. Paper search: K in [1,10], alpha in [0,5000].
% Numeric alpha must be positive in this implementation; log-search [10,5000].
cfg.searchK = [2, 8];                 % faster stage; change to [1,10] for full sweep
cfg.searchAlpha = [10, 5000];
cfg.fvim.population = 6;             % exploratory; increase to 12-20 for final runs
cfg.fvim.iterations = 4;             % exploratory; increase to 15-30 for final runs
cfg.mvmd.searchIterations = 65;
cfg.mvmd.finalIterations = 180;
cfg.mvmd.tolerance = 1e-6;
cfg.mvmd.tau = 0.025;

% Operationalized thresholds: NOT reported as unique values in the paper.
cfg.class.csKeep = 0.22;
cfg.class.sampKeep = 1.5;
cfg.class.csMixed = 0.08;
cfg.class.sampMixed = 2.2;
cfg.class.maxSampPoints = 320;
cfg.class.sampTolerance = 0.20;
cfg.wavelet.levels = 4;
cfg.wavelet.thresholdScale = 0.6;   % self-contained Haar-wavelet threshold

% Proposed (ours) physical constraints; NOT in the cited paper.
cfg.physics.ampTol = 0.01;           % 1% amplitude change
cfg.physics.phaseTolDeg = 0.20;      % degrees of two-channel relative phase change
cfg.physics.protectHz = 7;           % protects candidate modes near f0, not a hard signal cutoff
cfg.physics.minCarrierShare = 0.12;
cfg.physics.weights = [0.45, 0.35, 0.20]; % amplitude, phase, band energy
cfg.physics.bandHz = 12;            % sideband diagnostic span; not a universal cutoff

% Add genuine files, or use synthetic files generated at E00 for cross-frequency trial.
cfg.multiFiles = { ...
    fullfile(root, 'data', 'synthetic_dualgyro_3Hz.txt'), ...
    fullfile(root, 'data', 'synthetic_dualgyro_5Hz.txt'), ...
    fullfile(root, 'data', 'synthetic_dualgyro_10Hz.txt'), ...
    fullfile(root, 'data', 'synthetic_dualgyro_19Hz.txt')};
end
