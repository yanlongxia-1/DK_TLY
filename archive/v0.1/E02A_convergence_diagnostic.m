%% E02A: MVMD convergence and reconstruction diagnostic (non-destructive)
% Put this file next to gyro_config.m and gyro_mvmd.m.
% Requires E02/record.mat generated for the same cfg.file.
% No parameter fitting or threshold editing is performed here.
clear; clc; close all;
cfg = gyro_config();
[~,stem] = fileparts(cfg.file);
inputFile = fullfile(cfg.outputRoot,stem,'E02','record.mat');
if ~isfile(inputFile)
    error('E02 record.mat missing: %s\nSet cfg.file to your E02 source TXT and run E02 first.',inputFile);
end
S = load(inputFile,'record');
rec = S.record;
x = double(rec.raw);
K = double(rec.k);
alpha = double(rec.alpha);
fs = cfg.fs;
f0 = cfg.rollHz;
if size(x,2) ~= 2
    error('Validation signal must contain exactly two channels.');
end
fprintf('E02A: K=%d, alpha=%.8g, %d samples, fs=%.3f, rollHz=%.4f\n',K,alpha,size(x,1),fs,f0);

% First four experiments hold tau constant and only change the iteration budget.
% Next two test different dual-ascent step values, to diagnose slow/stalled convergence.
maxIterations = [180 360 720 1440 720 720];
tauValues = [cfg.mvmd.tau cfg.mvmd.tau cfg.mvmd.tau cfg.mvmd.tau 0.10 0.25];
nTests = numel(maxIterations);
T = table();
outDir = fullfile(cfg.outputRoot,stem,'E02A');
if ~isfolder(outDir), mkdir(outDir); end
allModes = cell(nTests,1);
allClasses = cell(nTests,1);

for j = 1:nTests
    imax = maxIterations(j);
    tau = tauValues(j);
    fprintf('\n[%d/%d] maxIter=%d, tau=%.4g ...\n',j,nTests,imax,tau);
    tt = tic;
    [modes,centerHz,info] = gyro_mvmd(x,fs,K,alpha,imax,cfg.mvmd.tolerance,tau);
    elapsed = toc(tt);
    [y,cl] = gyro_reconstruct(modes,centerHz,x,cfg,'paper');
    metrics = gyro_metrics(x,y,fs,f0,'E02A');
    xAC = x-mean(x,1);
    reconstructedAC = sum(modes,3);
    residualAC = xAC-reconstructedAC;
    relAC = norm(residualAC(:))/max(norm(xAC(:)),eps);
    fullRec = reconstructedAC+mean(x,1);
    relFull = norm(fullRec(:)-x(:))/max(norm(x(:)),eps);
    clValues = cl.class_0drop_1mix_2keep;
    row = table(imax,tau,info.iterations,double(info.converged),info.delta,...
        relAC,relFull,info.reconRelative,elapsed,...
        sum(clValues==0),sum(clValues==1),sum(clValues==2),...
        metrics.carrier_amp_keep_pct(1),metrics.carrier_amp_keep_pct(2),...
        metrics.phase_shift_deg(1),metrics.phase_shift_deg(2),...
        metrics.relative_phase_change_deg(1),...
        metrics.snr_proxy_gain_dB(1),metrics.snr_proxy_gain_dB(2),...
        metrics.mean_offset(1),metrics.mean_offset(2),...
        'VariableNames',{'requested_iterations','tau','actual_iterations','converged','final_delta',...
        'AC_reconstruction_relative','full_signal_reconstruction_relative','old_reconstruction_relative','runtime_seconds',...
        'n_discarded','n_mixed','n_kept','amp_keep_a_pct','amp_keep_s_pct',...
        'phase_shift_a_deg','phase_shift_s_deg','relative_phase_change_deg',...
        'snr_proxy_gain_a_dB','snr_proxy_gain_s_dB','mean_offset_a','mean_offset_s'});
    if j==1, T=row; else, T=[T;row]; end %#ok<AGROW>
    allModes{j} = centerHz(:);
    allClasses{j} = clValues(:);
    fprintf('  iter=%d, converged=%d, delta=%.3g, AC residual=%.4f%%, full residual=%.4f%%\n',...
        info.iterations,info.converged,info.delta,100*relAC,100*relFull);
    fprintf('  classification: discarded=%d, mixed=%d, kept=%d; phase delta=%.5f deg\n',...
        row.n_discarded,row.n_mixed,row.n_kept,row.relative_phase_change_deg);
end

writetable(T,fullfile(outDir,'E02A_convergence.csv'));
modeRows = table();
for j=1:nTests
    kcol=(1:numel(allModes{j}))';
    Srow=table(repmat(j,numel(kcol),1),kcol,allModes{j},allClasses{j},...
        'VariableNames',{'test_index','mode','center_Hz','class_0drop_1mix_2keep'});
    modeRows=[modeRows;Srow]; %#ok<AGROW>
end
writetable(modeRows,fullfile(outDir,'E02A_modes.csv'));

fig=figure('Visible','off','Position',[50 50 960 770]);
subplot(3,1,1);
plot(1:nTests,T.AC_reconstruction_relative*100,'-o','LineWidth',1.4); grid on;
ylabel('AC recon residual (%)'); xticks(1:nTests);title('Same original E02 validation data; K and alpha fixed');
subplot(3,1,2);
semilogy(1:nTests,max(T.final_delta,realmin),'-o','LineWidth',1.4);grid on;
ylabel('Final iteration relative change'); xticks(1:nTests);
subplot(3,1,3);
plot(1:nTests,T.snr_proxy_gain_a_dB,'-o',1:nTests,T.snr_proxy_gain_s_dB,'-s','LineWidth',1.4);
ylabel('SNR proxy gain (dB)');xlabel('Test index (see E02A_convergence.csv)');
legend('a','s','Location','best');grid on;xticks(1:nTests);
saveas(fig,fullfile(outDir,'E02A_diagnostic.png'));close(fig);

fid=fopen(fullfile(outDir,'E02A_README.txt'),'w');
fprintf(fid,['E02A compares MVMD convergence on the SAME E02 validation signal.\n',...
    'K and alpha remain fixed, no retuning.\n',...
    'AC_reconstruction_relative = norm(x-mean(x)-sum(IMFs))/norm(x-mean(x)).\n',...
    'This matters more than residual normalized by x including large sensor offsets.\n',...
    'E02A cannot establish true SNR or real angular-rate accuracy.\n',...
    'Do not choose optimum based on validation SNR or phase results; use diagnostic only.\n',...
    'Send E02A_convergence.csv, E02A_modes.csv, E02A_diagnostic.png.\n']);
fclose(fid);
fprintf('\nE02A finished. Outputs: %s\n',outDir);
disp(T(:,{'requested_iterations','tau','converged','final_delta','AC_reconstruction_relative',...
    'n_discarded','n_mixed','n_kept','relative_phase_change_deg'}));
