%% E03A - Fixed MVMD decomposition: check what physics-based IMF protection REALLY changes.
% Usage: Copy this file into the v0.3 MATLAB project root. Set cfg.file
% to the same TXT used for your successful E02. Run this script only.
%
% Crucial safeguards:
% 1) Reuses E02 record.mat; runs ZERO FVIM searches and ZERO MVMD decompositions.
% 2) Holds K, alpha, modes, validation indices, and all classification
%    thresholds fixed. Only switches the IMF reconstruction strategy.
% 3) Does NOT fit/choose new parameters on held-out data.
% 4) If no IMF is protected, reports exactly that; it does not manufacture
%    an improvement and cannot replace future real-truth validation.
clear; clc; close all;
startTimer=tic;
cfgCurrent = gyro_config();
[~, stem] = fileparts(cfgCurrent.file);
e02Path = fullfile(cfgCurrent.outputRoot, stem, 'E02', 'record.mat');
if ~isfile(e02Path)
    error('E03A:MissingE02', 'Cannot find E02 record.mat: %s', e02Path);
end
S = load(e02Path, 'record');
r = S.record;
required = {'raw','filt','modes','centers','cfg','k','alpha','valIdx','strategy','filename'};
for j=1:numel(required)
    if ~isfield(r, required{j})
        error('E03A:MissingField','E02 record lacks field: %s',required{j});
    end
end
if ~strcmpi(char(r.strategy), 'paper')
    error('E03A:WrongMethod','The input record must be the E02 paper baseline.');
end
if ~strcmp(char(r.filename), char(cfgCurrent.file))
    error('E03A:FileMismatch','cfg.file must match the E02 record filename.');
end
if size(r.raw,2)~=2 || size(r.modes,1)~=size(r.raw,1) || ...
        size(r.modes,2)~=2 || size(r.modes,3)~=r.k
    error('E03A:WrongShape','E02 record signal and mode array shapes disagree.');
end
if isfield(r,'mvmdInfo') && isfield(r.mvmdInfo,'converged') && ~r.mvmdInfo.converged
    error('E03A:Nonconverged','Saved E02 MVMD did not converge.');
end
% IMPORTANT: use the EXACT saved E02 config to avoid secretly changing
% threshold criteria or frequency while comparing reconstruction policies.
cfg = r.cfg;
[baselineCheck, tableBase] = gyro_reconstruct(r.modes,r.centers,r.raw,cfg,'paper');
[physicsFilt, tablePhysics] = gyro_reconstruct(r.modes,r.centers,r.raw,cfg,'physics');
repeatErr=max(abs(baselineCheck(:)-r.filt(:)));
if repeatErr > 1e-9 * max(1,max(abs(r.filt(:))))
    error('E03A:ReproductionMismatch', ...
        'Reconstructed E02 differs from stored E02 (max abs=%g). Do not compare.',repeatErr);
end
baselineClass=tableBase.class_0drop_1mix_2keep;
physicsClass=tablePhysics.class_0drop_1mix_2keep;
modeReport=tablePhysics;
modeReport.E02_class=baselineClass;
modeReport.E03A_class=physicsClass;
modeReport.class_changed=double(baselineClass~=physicsClass);
modeReport = movevars(modeReport,{'E02_class','E03A_class','class_changed'},'After','mode');
T02=gyro_metrics(r.raw,r.filt,cfg.fs,cfg.rollHz,'E02_frozen');
T03A=gyro_metrics(r.raw,physicsFilt,cfg.fs,cfg.rollHz,'E03A_fixed_MVMD_physics_selector');
allMetrics=[T02;T03A];
numChanged=sum(baselineClass~=physicsClass);
numProtected=sum(tablePhysics.physical_protection_applied>0);
maxAbsDelta=max(abs(physicsFilt(:)-r.filt(:)));
if isfield(r,'mvmdInfo') && isfield(r.mvmdInfo,'AC_reconRelative')
    acResidual=r.mvmdInfo.AC_reconRelative;
else
    acResidual=NaN;
end
summary=table(r.k,r.alpha,size(r.raw,1),numChanged,numProtected, ...
    maxAbsDelta,repeatErr,acResidual, ...
    'VariableNames',{'K_fixed','alpha_fixed','validation_samples','IMF_class_changed', ...
    'IMF_physics_protected','max_absolute_output_difference','E02_reproduction_maxdiff', ...
    'MVMD_AC_recon_relative'});
out=fullfile(cfgCurrent.outputRoot,stem,'E03A');
if ~isfolder(out),mkdir(out);end
writetable(summary,fullfile(out,'E03A_summary.csv'));
writetable(modeReport,fullfile(out,'E03A_mode_changes.csv'));
writetable(allMetrics,fullfile(out,'E03A_filter_metrics.csv'));
t=(r.valIdx(:)-1)/cfg.fs;
writematrix([t,r.raw,r.filt,physicsFilt], ...
    fullfile(out,'E03A_validation_compare.csv'));
fig=figure('Visible','off','Color','w');
labels={'a / column 4','s / column 7'};
for c=1:2
    subplot(2,1,c);
    plot(t,r.raw(:,c),'Color',[0.75 0.75 0.75]);hold on;
    plot(t,r.filt(:,c),'b','LineWidth',1);
    plot(t,physicsFilt(:,c),'r--','LineWidth',1);
    legend('Original','E02 baseline','E03A fixed-modes physics','Location','best');
    title(labels{c});ylabel('Sensor units');grid on;
end
xlabel('Time (s)');
saveas(fig,fullfile(out,'E03A_comparison.png'));close(fig);
fid=fopen(fullfile(out,'E03A_readme.txt'),'w');
fprintf(fid,'E03A is a FIXED-MVMD, FIXED-PARAMETER selector-only ablation.\n');
fprintf(fid,'No FVIM or MVMD search/decomposition was executed in E03A.\n');
fprintf(fid,'Only saved E02 validation modes were reconstructed twice.\n');
fprintf(fid,'This is a DIAGNOSTIC on an already viewed validation interval, not a fresh independent efficacy test.\n');
fprintf(fid,'No measured angular-rate truth is available; SNR values are fit proxies.\n');
fprintf(fid,'DO NOT choose classification thresholds using these validation results.\n');
fprintf(fid,'Zero protected IMFs means the current physical guard had no effect on this signal.\n');
fprintf(fid,'K=%d, alpha=%.8f, changed_classes=%d, physical_protection=%d, max_output_difference=%.12g\n', ...
    r.k,r.alpha,numChanged,numProtected,maxAbsDelta);
fclose(fid);
fprintf('\nE03A completed in %.2f s; K=%d, alpha=%.6f\n',toc(startTimer),r.k,r.alpha);
fprintf('IMF class changes=%d; physics protections=%d; max output change=%.8g\n', ...
    numChanged,numProtected,maxAbsDelta);
if numChanged==0
    fprintf('FINDING: Physical IMF guarding does not change reconstruction here.\n');
    fprintf('Next scientific step: dynamic sideband / ground-truth experiment, not threshold tuning on validation.\n');
end
fprintf('Outputs: %s\n',out);
