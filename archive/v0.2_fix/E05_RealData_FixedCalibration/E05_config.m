function cfg = E05_config()
% E05: REAL W3-10 / W3-20, fixed per-condition bias and scale.
% Set the 8 NaN calibration fields below to independently measured numbers.
% DO NOT estimate them from the validation interval.
root = fileparts(mfilename('fullpath'));
cfg.fs = 1000;
cfg.outputRoot = fullfile(root,'results');
cfg.trainFraction = [0.10,0.50];
cfg.validationFraction = [0.60,0.95];

C = struct('file',{},'fRoll',{},'phi0Deg',{},'rollSign',{},'bias',{},'scale',{});
C(1).file = fullfile(root,'data','W3-10.txt');
C(1).fRoll = 3.0; C(1).phi0Deg = 0; C(1).rollSign = +1;
C(1).bias = [NaN NaN];       % [b_a,b_s] in recorded sensor-output units
C(1).scale = [NaN NaN];      % [K_a,K_s] sensor-output units / (deg/s)
C(2).file = fullfile(root,'data','W3-20.txt');
C(2).fRoll = 3.0; C(2).phi0Deg = 0; C(2).rollSign = +1;
C(2).bias = [NaN NaN];
C(2).scale = [NaN NaN];
cfg.cases = C;
% phi(t) = phi0 + rollSign*2*pi*fRoll*t; identical for all methods.
% phi0=0 and rollSign=+1 are explicit assumptions, NOT fitted results.

% FROZEN E02 baseline parameters originally selected with W3-10.
% Not reoptimised on validation or on W3-20.
cfg.mvmd.K = 8;
cfg.mvmd.alpha = 154.1579598620;
cfg.mvmd.maxIterations = 1200;
cfg.mvmd.tolerance = 1e-5;
cfg.mvmd.tau = 0;

% Fixed KF baseline; train-only R estimate, no validation-time tuning.
% This uses the same oscillator structure, not a byte-for-byte reproduction
% of any separately tuned older KF program.
cfg.kf.qOscOverR = 1e-6;
cfg.kf.qBiasOverR = 1e-8;

% FROZEN E02 paper-style operationalization (not original author source).
cfg.class.csKeep = 0.22;
cfg.class.sampKeep = 1.5;
cfg.class.csMixed = 0.08;
cfg.class.sampMixed = 1.8;
cfg.class.maxSampPoints = 320;
cfg.class.sampTolerance = 0.2;
cfg.wavelet.levels = 4;
cfg.wavelet.thresholdScale = 0.6;

% Reference attitudes from col1/2/3: currently disabled.
% To enable, physical gimbal axes/order, signs, sensor mounting, and
% meaning of Omega_x/Omega_y MUST have been independently confirmed.
% Do not switch this on based only on a low RMSE result.
cfg.reference.enabled = false;
cfg.reference.axes = {'UNVERIFIED','UNVERIFIED','UNVERIFIED'};
cfg.reference.signs = [NaN NaN NaN];
cfg.reference.sensorMount = eye(3);
cfg.reference.omegaXYFrame = 'UNVERIFIED'; % only 'body' accepted by current module
cfg.reference.rollBodyAxis = 'UNVERIFIED'; % only 'z' accepted by current module
end
