function E03E_phase_reference_and_calibration()
% E03E | Phase observability, reference requirement, and separate calibration errors
% Self-contained MATLAB experiment. No FVIM, MVMD, toolbox, or external data.
% IMPORTANT: All roll-rate and absolute-phase references below are SIMULATED.
% They are NOT claims about a real accelerometer or the 6th TXT column.
% This experiment does not measure real-sensor accuracy.
clearvars; close all; clc;

fs = 1000; N = 2048; t = (0:N-1)'/fs; T=t(end);
f0 = 3; df = 0.10; phi0 = 5*pi/180;
Ka0 = 0.0020; Ks0 = 0.0022;
b0 = [-0.186041,-0.072559];
Ox = 20 + 8*sin(2*pi*3*t);
Oy = 8 + 6*cos(2*pi*3*t);
roi = (101:N-100)';

% Controlled phase-error case: both an unknown 5-degree initial phase
% and a 3.00 -> 3.10 Hz linear ramp over T seconds.
ftrue = f0 + df*t/T;
phiTrue = 2*pi*(f0*t + 0.5*df*t.^2/T) + phi0;
a = Ka0*(Ox.*sin(phiTrue) - Oy.*cos(phiTrue)) + b0(1);
s = Ks0*(Ox.*cos(phiTrue) + Oy.*sin(phiTrue)) + b0(2);

% Check the exact inverse before any phase/reference experiments.
[tx,ty] = inverse_rates(a,s,phiTrue,b0,[Ka0,Ks0]);
checkError = max([abs(tx-Ox); abs(ty-Oy)]);
assert(checkError < 1e-8, 'E03E: inverse-model self-check failed.');

outDir = fullfile(pwd,'E03E_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
numSeeds = 20;
seeds = 20271010+(1:numSeeds);
methods = {'true_phi_oracle','fixed_3Hz_zero_phase', ...
    'external_roll_rate_only','external_rate_known_initial_phase', ...
    'rate_plus_external_absolute_phase'};

% Hypothetical independent roll-rate source: noise std 0.02 Hz,
% constant bias +0.005 Hz. The 5-degree initial phase is NOT observable
% from the rate source alone.
rateNoiseHz = 0.02; rateBiasHz = 0.005;
% Hypothetical independent absolute-phase reference: every 0.25 seconds,
% noise std 2 degrees. This is NOT an accelerometer performance assumption.
anchorStep = round(0.25*fs);
anchorNoiseDeg = 2;
anchorGain = 0.7;
anchors = 1:anchorStep:N;

phaseRows = cell(numSeeds*numel(methods),9);
r = 0;
for j=1:numSeeds
    rng(seeds(j),'twister');
    fr = ftrue + rateBiasHz + rateNoiseHz*randn(N,1);
    % cumtrapz avoids a spurious initial integration step.
    phiRate = 2*pi*cumtrapz(t,fr);
    phiFixed = 2*pi*f0*t;
    phiRateKnownInit = phiRate + phi0; % informative oracle initial anchor
    phiFusion = zeros(N,1);
    isAnchor = false(N,1); isAnchor(anchors) = true;
    for n=1:N
        if n>1
            phiFusion(n) = phiFusion(n-1) + 2*pi*0.5*(fr(n-1)+fr(n))/fs;
        end
        if isAnchor(n)
            absoluteObservation = phiTrue(n) + (anchorNoiseDeg*pi/180)*randn();
            e = angular_difference(absoluteObservation,phiFusion(n));
            phiFusion(n) = phiFusion(n) + anchorGain*e;
        end
    end
    phaseEstimates = {phiTrue,phiFixed,phiRate,phiRateKnownInit,phiFusion};
    for k=1:numel(methods)
        ph=phaseEstimates{k};
        [xx,yy] = inverse_rates(a,s,ph,b0,[Ka0,Ks0]);
        pe = angular_difference(ph(roi),phiTrue(roi))*180/pi;
        pend = angular_difference(ph(end),phiTrue(end))*180/pi;
        ex = sqrt(mean((xx(roi)-Ox(roi)).^2));
        ey = sqrt(mean((yy(roi)-Oy(roi)).^2));
        em = sqrt(mean((hypot(xx(roi),yy(roi))-hypot(Ox(roi),Oy(roi))).^2));
        r=r+1;
        phaseRows(r,:)={'offset5_plus_ramp0p10',seeds(j),methods{k}, ...
            sqrt(mean(pe.^2)),pend,ex,ey,hypot(ex,ey),em};
    end
end
phaseTable=cell2table(phaseRows,'VariableNames',{'scenario','seed','method', ...
    'phase_RMSE_deg','phase_end_error_deg','OmegaX_RMSE_deg_s', ...
    'OmegaY_RMSE_deg_s','OmegaVector_RMSE_deg_s','OmegaMagnitude_RMSE_deg_s'});
writetable(phaseTable,fullfile(outDir,'E03E_phase_seed_results.csv'));

summ=cell(numel(methods),7);
for k=1:numel(methods)
    sub=phaseTable(strcmp(phaseTable.method,methods{k}),:);
    summ(k,:)={methods{k},mean(sub.phase_RMSE_deg), ...
        mean(abs(sub.phase_end_error_deg)),mean(sub.OmegaX_RMSE_deg_s), ...
        mean(sub.OmegaY_RMSE_deg_s),mean(sub.OmegaVector_RMSE_deg_s), ...
        mean(sub.OmegaMagnitude_RMSE_deg_s)};
end
summaryTable=cell2table(summ,'VariableNames',{'method','mean_phase_RMSE_deg', ...
    'mean_abs_phase_end_error_deg','mean_OmegaX_RMSE_deg_s', ...
    'mean_OmegaY_RMSE_deg_s','mean_OmegaVector_RMSE_deg_s', ...
    'mean_OmegaMagnitude_RMSE_deg_s'});
writetable(summaryTable,fullfile(outDir,'E03E_phase_summary.csv'));

% Independent zero-noise calibration experiment: vary bias and gains in
% separate conditions so E03D's coupled case is not mistaken for attribution.
phiC = 2*pi*f0*t;
labels={'correct','bias_only','gain_only','bias_and_gain'};
calRows=cell(numel(labels),8);
for k=1:numel(labels)
    addBias = (k==2 || k==4);
    addGain = (k==3 || k==4);
    bTrue=b0; KTrue=[Ka0,Ks0];
    if addBias, bTrue=bTrue+[0.004,-0.003]; end
    if addGain, KTrue=KTrue.*[1.06,0.94]; end
    aC=KTrue(1)*(Ox.*sin(phiC)-Oy.*cos(phiC))+bTrue(1);
    sC=KTrue(2)*(Ox.*cos(phiC)+Oy.*sin(phiC))+bTrue(2);
    [xx,yy]=inverse_rates(aC,sC,phiC,b0,[Ka0,Ks0]);
    ex=sqrt(mean((xx(roi)-Ox(roi)).^2));
    ey=sqrt(mean((yy(roi)-Oy(roi)).^2));
    em=sqrt(mean((hypot(xx(roi),yy(roi))-hypot(Ox(roi),Oy(roi))).^2));
    calRows(k,:)={labels{k},bTrue(1)-b0(1),bTrue(2)-b0(2), ...
        100*(KTrue(1)/Ka0-1),100*(KTrue(2)/Ks0-1),ex,ey,em};
end
calTable=cell2table(calRows,'VariableNames',{'condition', ...
    'actual_bias_mismatch_a','actual_bias_mismatch_s', ...
    'actual_gain_mismatch_a_pct','actual_gain_mismatch_s_pct', ...
    'OmegaX_RMSE_deg_s','OmegaY_RMSE_deg_s','OmegaMagnitude_RMSE_deg_s'});
writetable(calTable,fullfile(outDir,'E03E_calibration_separate.csv'));

% Light-weight explanatory figure. No heavy decomposition involved.
fig=figure('Visible','off','Color','w','Position',[100 100 900 650]);
subplot(2,1,1);
bar([summaryTable.mean_phase_RMSE_deg]);
set(gca,'XTick',1:numel(methods),'XTickLabel',methods);
xtickangle(18); ylabel('Phase RMSE (deg)'); grid on;
title('Synthetic external phase-reference requirements');
subplot(2,1,2);
bar([summaryTable.mean_OmegaX_RMSE_deg_s,summaryTable.mean_OmegaY_RMSE_deg_s]);
set(gca,'XTick',1:numel(methods),'XTickLabel',methods);
xtickangle(18); ylabel('Signed rate RMSE (deg/s)');
legend('OmegaX','OmegaY','Location','best'); grid on;
try, exportgraphics(fig,fullfile(outDir,'E03E_phase_methods.png'),'Resolution',140);
catch, saveas(fig,fullfile(outDir,'E03E_phase_methods.png')); end
close(fig);

fid=fopen(fullfile(outDir,'E03E_notes.txt'),'wt');
fprintf(fid,'E03E: analytic inverse max absolute error = %.4g deg/s\n',checkError);
fprintf(fid,'Pure source-isolation experiments: CLEAN sensor signals, NO simulated ADC noise, NO MVMD or FVIM.\n');
fprintf(fid,'Scenario: unknown initial phase = +5 deg, roll frequency ramps 3.00 -> 3.10 Hz over 2.047 seconds.\n');
fprintf(fid,'External roll-rate observation: bias +0.005 Hz, Gaussian noise std 0.02 Hz (SIMULATED).\n');
fprintf(fid,'Absolute phase observation: every 0.25s, noise std 2deg, fusion gain 0.7 (SIMULATED).\n');
fprintf(fid,'Known-initial-phase method is an ORACLE calibration reference, not a directly available method.\n');
fprintf(fid,'Bias-only and gain-only are isolated in the noiseless second experiment.\n');
fprintf(fid,'Magnitude invariance to phase rotation is expected; check signed rates and vector RMSE.\n');
fprintf(fid,'Independent phase reference is hypothetical; col6 accelerometer or turntable column3 NOT validated.\n');
fprintf(fid,'This is a requirements study, NOT evidence that real IMU/accelerometer yields phase reference.\n');
fprintf(fid,'Generated in: %s\n',outDir);
fclose(fid);

disp(summaryTable); disp(calTable);
fprintf('\nE03E complete. Output folder: %s\n',outDir);
end

function [ox,oy]=inverse_rates(a,s,phi,b,K)
a1=(a-b(1))/K(1); s1=(s-b(2))/K(2);
ox=s1.*cos(phi)+a1.*sin(phi);
oy=s1.*sin(phi)-a1.*cos(phi);
end

function d=angular_difference(a,b)
d=atan2(sin(a-b),cos(a-b));
end
