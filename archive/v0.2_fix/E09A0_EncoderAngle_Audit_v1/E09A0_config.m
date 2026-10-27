function cfg = E09A0_config()
% E09A0: 不重新滤波；仅用E08-K固定输出与原始第1列角度，验证积分角度。
root = fileparts(mfilename('fullpath'));
cfg.dataRoot = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验';
cfg.encoderTxtName = 'W9-60.TXT';
cfg.rateCsv = fullfile(root, '..', 'E08K_AutoFVIM_Physics_Rate_v1', ...
    'results_E08K_9Hz','E08K_series_W9_60_TXT.csv');
cfg.outputDir = fullfile(root, 'results_E09A0');
cfg.fs = 1000;
cfg.sampleTimeToleranceSec = 1e-7;
cfg.referenceAngleCol = 1; % 预先确认此列为转台外框角度，单位deg
cfg.makeFigure = true;
end
