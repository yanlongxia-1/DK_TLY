%% E02C: train-only threshold sensitivity, DO NOT tune from validation metrics
% Requires corrected E02 record.mat (converged) to fix the K and alpha only.
clear;clc;close all;
cfg=gyro_config();
[~,stem]=fileparts(cfg.file);
recFile=fullfile(cfg.outputRoot,stem,'E02','record.mat');
if ~isfile(recFile),error('Please run corrected v0.2 E02 first: %s',recFile);end
S=load(recFile,'record');r=S.record;
if ~isfield(r,'version') || ~strcmp(r.version,'v0.2')
    error('Old v0.1 record detected. Run corrected E02 first, do not mix versions.');
end
d=gyro_read(cfg.file,cfg);
x=d.y(r.trainIdx,:);
[modes,cent,info]=gyro_mvmd(x,cfg.fs,r.k,r.alpha,...
    cfg.mvmd.finalIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
if ~info.converged,error('Training MVMD failed convergence');end
out=fullfile(cfg.outputRoot,stem,'E02C');if ~isfolder(out),mkdir(out);end
csOptions=[0.06,0.08,0.10,0.12];
spOptions=[1.5,1.8,2.2];
rows={};
for ic=1:numel(csOptions)
  for is=1:numel(spOptions)
    cc=cfg;cc.class.csMixed=csOptions(ic);cc.class.sampMixed=spOptions(is);
    [y,dt]=gyro_reconstruct(modes,cent,x,cc,'paper');
    met=gyro_metrics(x,y,cfg.fs,cfg.rollHz,'E02C_train');
    rows(end+1,:)={cc.class.csMixed,cc.class.sampMixed,...
       sum(dt.class_0drop_1mix_2keep==0),...
       sum(dt.class_0drop_1mix_2keep==1),...
       sum(dt.class_0drop_1mix_2keep==2),...
       max(abs(met.carrier_amp_keep_pct-100)),...
       abs(met.relative_phase_change_deg(1)),...
       met.snr_proxy_gain_dB(1),met.snr_proxy_gain_dB(2)}; %#ok<AGROW>
  end
end
T=cell2table(rows,'VariableNames',{'cs_mixed','sampen_mixed',...
    'n_dropped','n_mixed','n_kept','max_carrier_amp_deviation_pct',...
    'relative_phase_change_deg','snr_proxy_gain_a_dB','snr_proxy_gain_s_dB'});
writetable(T,fullfile(out,'E02C_train_thresholds.csv'));
writetable(gyro_reconstruct_details(modes,cent,x,cfg),...
    fullfile(out,'E02C_fixed_rule_modes.csv'));
fid=fopen(fullfile(out,'E02C_notes.txt'),'w');
fprintf(fid,'Only train interval used. Changing thresholds based on validation data is forbidden.\n');
fprintf(fid,'CS and SampEn thresholds are ENGINEERING operationalization; paper has no unique universal threshold.\n');
fprintf(fid,'Training SNR is a single-tone PROXY, not proof of true noise removal.\n');
fprintf(fid,'Any thresholds changed after looking at the old validation outcomes need testing on a NEW independent holdout.\n');
fclose(fid);
disp(T);fprintf('E02C saved: %s\n',out);
function tab=gyro_reconstruct_details(m,cent,x,cfg)
 [~,tab]=gyro_reconstruct(m,cent,x,cfg,'paper');
end
