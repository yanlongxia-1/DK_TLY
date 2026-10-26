
clear; clc; close all;

cfg = E07B_config();

% 只读取训练数据，不访问Z3-2
D = load(cfg.inputTrain, 'X', 't');

alphaList = [292.97310, 306.12870];
K = 8;

for i = 1:numel(alphaList)
    alpha = alphaList(i);

    fprintf('\n===== K=%d, alpha=%.5f =====\n', K, alpha);

    [modes, centerHz, info] = gyro_mvmd( ...
        D.X, cfg.fs, K, alpha, ...
        cfg.mvmd.maxIterations, ...
        cfg.mvmd.tolerance, cfg.mvmd.tau);

    fprintf('收敛状态: %d\n', info.converged);
    fprintf('迭代次数: %d\n', info.iterations);
    fprintf('AC重构误差: %.6f\n', info.AC_reconRelative);

    [Y, details] = gyro_reconstruct( ...
        modes, centerHz, D.X, cfg, 'paper');

    metrics = E07B_candidate_metrics( ...
        D.X, Y, D.t, ...
        cfg.trainComponentsHz, ...
        cfg.trainSidebandHz, cfg);

    disp(details(:, { ...
        'mode', 'center_Hz', ...
        'cosine_abs', 'sampEn', ...
        'energy_share', ...
        'class_0drop_1mix_2keep'}));

    fprintf('a通道降噪代理: %.4f dB\n', ...
        metrics.nonTargetSuppressionDb(1));
    fprintf('s通道降噪代理: %.4f dB\n', ...
        metrics.nonTargetSuppressionDb(2));

    name = sprintf('E07B_IMF_alpha%d.csv', round(alpha));
    writetable(details, fullfile(cfg.outputDir, name));
end

disp('诊断完成：已导出两组IMF明细。');
