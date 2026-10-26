%% E03B: ground-truth stress test, FIXED E02 K/alpha, no FVIM search.
% Copy to v0.3 project root. Change only gyro_config.m cfg.file to your E02 TXT.
% Stresses stable carrier, dynamic transverse angular rates, and roll-speed drift.
% Outputs true sensor errors and TRUE signed transverse rate RMSE for raw/E02/physics.
clear; clc; close all;
cfgCurrent=gyro_config();
[~,stem]=fileparts(cfgCurrent.file);
p=fullfile(cfgCurrent.outputRoot,stem,'E02','record.mat');
if ~isfile(p),error('E03B:NeedE02','E02 record.mat missing: %s',p);end
Q=load(p,'record');R=Q.record;
if ~isfield(R,'mvmdInfo') || ~R.mvmdInfo.converged
    error('E03B:NeedConvergedE02','E02 saved solver was not converged.');
end
% Freeze all choices from E02, no future test leakage.
cfg=R.cfg;
N=2048;fs=cfg.fs;t=(0:N-1)'/fs;
Ka=0.0020;Ks=0.0022;
b=[-0.186041,-0.072559];
sigma=0.012;
runSeed=20261025+301;
rng(runSeed,'twister');
caseNames={'steady_3Hz','varying_OmegaXY','varying_roll_and_OmegaXY'};
out=fullfile(cfgCurrent.outputRoot,stem,'E03B');
if ~isfolder(out),mkdir(out);end
rows={};
for scenario=1:3
    switch scenario
        case 1
            phi=2*pi*3*t;ox=20+zeros(N,1);oy=8+zeros(N,1);
        case 2
            phi=2*pi*3*t;
            ox=20+8*sin(2*pi*5*t);oy=8+6*cos(2*pi*6*t);
        case 3
            % instantaneous roll increases smoothly from ~3Hz to ~4.5Hz
            phi=2*pi*(3*t+0.5*0.75*t.^2)+0.07*sin(2*pi*0.45*t);
            ox=20+8*sin(2*pi*5*t);oy=8+6*cos(2*pi*6*t);
    end
    clean=[Ka*(ox.*sin(phi)-oy.*cos(phi)), ...
           Ks*(ox.*cos(phi)+oy.*sin(phi))];
    y=clean+b+sigma*randn(N,2);
    % Identical sparse glitches across filters (difficult-but-controlled).
    glitches=[402,1101,1719];
    y(glitches,1)=y(glitches,1)+[0.055;-0.047;0.05];
    y(glitches+9,2)=y(glitches+9,2)+[-0.052;0.041;-0.056];
    tic;
    [modes,cent,info]=gyro_mvmd(y,fs,R.k,R.alpha, ...
        min(1200,cfg.mvmd.finalIterations),cfg.mvmd.tolerance,cfg.mvmd.tau);
    runtime=toc;
    if info.converged
        [yp,Tp]=gyro_reconstruct(modes,cent,y,cfg,'paper');
        [yh,Th]=gyro_reconstruct(modes,cent,y,cfg,'physics');
        nProtect=sum(Th.physical_protection_applied);
        classChanged=sum(Tp.class_0drop_1mix_2keep~=Th.class_0drop_1mix_2keep);
    else
        yp=nan(size(y));yh=nan(size(y));nProtect=NaN;classChanged=NaN;
        warning('E03B:NonConverged','Scenario %s nonconverged, marked INVALID.',caseNames{scenario});
    end
    yy={y,yp,yh};labels={'raw','E02_paper_selector','E03_physics_selector'};
    for j=1:3
        v=yy{j};
        err=(v-b)-clean;
        sensorRMSE=sqrt(mean(err.^2,1));
        aa=(v(:,1)-b(1))/Ka; ss=(v(:,2)-b(2))/Ks;
        oxhat=ss.*cos(phi)+aa.*sin(phi);
        oyhat=ss.*sin(phi)-aa.*cos(phi);
        rateErr=sqrt(mean(([oxhat,oyhat]-[ox,oy]).^2,1));
        trueMagnitude=hypot(ox,oy);estimatedMagnitude=hypot(oxhat,oyhat);
        magnitudeRMSE=sqrt(mean((estimatedMagnitude-trueMagnitude).^2));
        rows(end+1,:)={caseNames{scenario},labels{j},double(info.converged), ...
            info.iterations,info.AC_reconRelative,runtime,R.k,R.alpha, ...
            sensorRMSE(1),sensorRMSE(2),rateErr(1),rateErr(2), ...
            magnitudeRMSE,nProtect,classChanged}; %#ok<AGROW>
    end
    T=cell2table(rows,'VariableNames',{'scenario','method','converged', ...
        'iterations','AC_recon_relative','MVMD_seconds','K_fixed','alpha_fixed', ...
        'TRUE_sensor_RMSE_a','TRUE_sensor_RMSE_s', ...
        'TRUE_OmegaX_RMSE_deg_s','TRUE_OmegaY_RMSE_deg_s', ...
        'TRUE_OmegaMagnitude_RMSE_deg_s','nProtectedIMFs','nChangedIMFClasses'});
    writetable(T,fullfile(out,'E03B_true_RMSE.csv'));
    if info.converged
        fprintf('%s: MVMD iterations=%d runtime=%.2fs; protection=%d, class changes=%d\n', ...
            caseNames{scenario},info.iterations,runtime,nProtect,classChanged);
    end
    if scenario>=2
        fig=figure('Visible','off','Color','w');
        subplot(2,1,1);plot(t,ox,'k','LineWidth',1.2);hold on;
        for j=1:3
            v=yy{j};aa=(v(:,1)-b(1))/Ka;ss=(v(:,2)-b(2))/Ks;
            plot(t,ss.*cos(phi)+aa.*sin(phi),'LineWidth',0.8);
        end
        title([caseNames{scenario},'  |  signed OmegaX']);grid on;
        legend('Truth','Raw','E02','E03 physical','Location','best');ylabel('deg/s');
        subplot(2,1,2);plot(t,oy,'k','LineWidth',1.2);hold on;
        for j=1:3
            v=yy{j};aa=(v(:,1)-b(1))/Ka;ss=(v(:,2)-b(2))/Ks;
            plot(t,ss.*sin(phi)-aa.*cos(phi),'LineWidth',0.8);
        end
        title('signed OmegaY');grid on;ylabel('deg/s');xlabel('Time (s)');
        saveas(fig,fullfile(out,[caseNames{scenario},'_signed_rates.png']));close(fig);
    end
end
fid=fopen(fullfile(out,'E03B_notes.txt'),'w');
fprintf(fid,'Synthetic reference is exact clean sensor output and signed OmegaX/OmegaY.\n');
fprintf(fid,'Known phi, bias, and sensitivities are used ONLY for simulation evaluation.\n');
fprintf(fid,'Both filters use the SAME frozen E02 K/alpha and SAME MVMD modes.\n');
fprintf(fid,'This isolates IMF selection; it DOES NOT evaluate physical FVIM parameter optimization.\n');
fprintf(fid,'No parameters may be tuned from E03B and then assessed on the same synthetic realization.\n');
fprintf(fid,'K=%d alpha=%.8g sigma=%.5g rngseed=%d\n',R.k,R.alpha,sigma,runSeed);
fclose(fid);
fprintf('E03B finished: %s\n',out);
disp(T);
