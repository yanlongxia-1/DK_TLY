function cfg = E05_calibrated_config
% E05 calibrated W bench validation: all fixed model coefficients are TRAIN-only.
root=fileparts(mfilename('fullpath'));
% Only field you MUST change: folder recursively containing E*.TXT / W*.TXT.
cfg.dataRoot='E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.calibrationCsv=fullfile(root,'calibration_final.csv');
cfg.outputRoot=fullfile(root,'results_E05_calibrated');
% Start with ONE case, then extend to e.g. {'W3-20.TXT','W3-60.TXT','W5-20.TXT'}
% W20/W60 were HELD OUT from the CALIBRATION regression, but prior filter
% inspection of W3-20 means that W3-20 is not a new blind FILTER test.
cfg.testFiles={'W3-20.TXT'};
cfg.fs=1000;
cfg.trainFraction=[0.10 0.50];
cfg.validationFraction=[0.60 0.95];
cfg.settlingSeconds=0.5;
cfg.endTrimSeconds=0.2;
cfg.minimumUsableSeconds=2.0;
cfg.huberTune=1.345;
cfg.huberMaxIterations=50;
cfg.huberConvergenceTolerance=1e-10;
cfg.outlierMadThreshold=6;
cfg.confidenceZ=1.96;
% FROZEN: original lab's MVMD parameters from W3-10. No FVIM.
cfg.mvmd.K=8;
cfg.mvmd.alpha=154.1579598620;
cfg.mvmd.maxIterations=1200;
cfg.mvmd.tolerance=1e-5;
cfg.mvmd.tau=0;
cfg.kf.qOscOverR=1e-6;
cfg.kf.qBiasOverR=1e-8;
cfg.class.csKeep=0.22;
cfg.class.sampKeep=1.5;
cfg.class.csMixed=0.08;
cfg.class.sampMixed=1.8;
cfg.class.maxSampPoints=320;
cfg.class.sampTolerance=0.2;
cfg.wavelet.levels=4;
cfg.wavelet.thresholdScale=0.6;
end
