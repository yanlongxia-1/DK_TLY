%% E05: SAME-real-data raw / fixed oscillator KF / fixed MVMD comparison.
% Running this script NEVER executes FVIM; no automatic calibration tuning.
clear; clc; close all;
cfg = E05_config();
if ~exist(cfg.outputRoot,'dir'),mkdir(cfg.outputRoot);end
for caseID=1:numel(cfg.cases)
    cs=cfg.cases(caseID);
    fprintf('\n========= E05 %s =========\n',cs.file);
    if ~isfile(cs.file)
        warning('Cannot find %s. Copy TXT into data/ or edit E05_config.m. Skipped.',cs.file);
        continue
    end
    try
        runOne(cs,cfg);
    catch ME
        [~,stem]=fileparts(cs.file);
        out=fullfile(cfg.outputRoot,stem);
        if ~isfolder(out),mkdir(out);end
        fid=fopen(fullfile(out,'E05_ERROR.txt'),'w');
        if fid>0
            fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));
            fclose(fid);
        end
        warning('E05 %s failed: %s',stem,ME.message);
    end
end
fprintf('\nE05 done. Look in results/ for individual-case reports.\n');

function runOne(cs,cfg)
    [~,stem]=fileparts(cs.file);
    out=fullfile(cfg.outputRoot,stem);
    if ~isfolder(out),mkdir(out);end
    M=readmatrix(cs.file,'FileType','text');
    if size(M,2)<7,error('Expected >=7 columns; got %d',size(M,2));end
    if any(~isfinite(M(:,[1,2,3,4,7])),'all')
        error('Invalid NaN/Inf in reference angles or gyro channels; do not silently delete samples.');
    end
    y=double(M(:,[4,7])); n=size(y,1);
    if n<1000,error('Too few samples for split: %d',n);end
    t=(0:n-1)'/cfg.fs;
    tr=clampIndices(n,cfg.trainFraction);
    va=clampIndices(n,cfg.validationFraction);
    if ~isempty(intersect(tr,va)),error('Train and validation overlap');end
    x=y(va,:);tv=t(va);
    fprintf('samples=%d; train=%d, validation=%d; frozen K=%d alpha=%.6f\n',...
        n,numel(tr),numel(va),cfg.mvmd.K,cfg.mvmd.alpha);
    % Baseline 1: raw measurements.
    % Baseline 2: fixed-Q oscillator KF; R and initialisation from TRAIN only.
    tic; kfAll=E05_KF_fixed(y,cfg.fs,cs.fRoll,tr,...
           cfg.kf.qOscOverR,cfg.kf.qBiasOverR);secKF=toc;
    yKF=kfAll(va,:);
    % Baseline 3: one frozen MVMD run, NO FVIM or refit.
    tic;
    [modes,centers,info]=gyro_mvmd(x,cfg.fs,cfg.mvmd.K,...
       cfg.mvmd.alpha,cfg.mvmd.maxIterations,...
       cfg.mvmd.tolerance,cfg.mvmd.tau);
    secMVMD=toc;
    fprintf('MVMD convergence=%d; iter=%d; AC residual=%.5f%%; time=%.2fs\n',...
       info.converged,info.iterations,100*info.AC_reconRelative,secMVMD);
    filtReady=logical(info.converged);
    if filtReady
        rcfg=cfg;rcfg.rollHz=cs.fRoll;
        [yMVMD,modeTable]=gyro_reconstruct(modes,centers,x,rcfg,'paper');
        writetable(modeTable,fullfile(out,'E05_modes.csv'));
    else
        yMVMD=nan(size(x));
        warning('MVMD did not converge. Will NOT report its signal/angle metrics.');
    end
    % Real gimbal angular POSITIONS, never gimbal rates.
    angle0=M(:,1:3);
    angleUnwrap=zeros(size(angle0));
    for j=1:3
        angleUnwrap(:,j)=rad2deg(unwrap(deg2rad(angle0(:,j))));
    end
    Tangle=table(t,angle0(:,1),angle0(:,2),angle0(:,3),...
      angleUnwrap(:,1),angleUnwrap(:,2),angleUnwrap(:,3),...
      'VariableNames',{'time_s','outer_deg','middle_deg','inner_deg',...
          'outer_unwrapped_deg','middle_unwrapped_deg','inner_unwrapped_deg'});
    writetable(Tangle,fullfile(out,'E05_gimbal_angles_position.csv'));
    % Same validation range and frequency-fit model for ALL methods.
    TT=gyro_metrics(x,x,cfg.fs,cs.fRoll,'raw');
    TT=[TT;gyro_metrics(x,yKF,cfg.fs,cs.fRoll,'KF_fixed')];
    if filtReady
        TT=[TT;gyro_metrics(x,yMVMD,cfg.fs,cs.fRoll,'FVIM_MVMD_frozen')];
    end
    writetable(TT,fullfile(out,'E05_signal_metrics.csv'));
    casestring={stem};
    S=table(casestring,cs.fRoll,cs.rollSign,cs.phi0Deg,...
      cfg.mvmd.K,cfg.mvmd.alpha,double(info.converged),...
      info.iterations,info.delta,info.AC_reconRelative,...
      secKF,secMVMD,numel(tr),numel(va),...
      'VariableNames',{'case_name','roll_frequency_Hz','roll_sign_assumed',...
       'phi0_deg_assumed','K_frozen','alpha_frozen',...
       'MVMD_converged','MVMD_iterations','MVMD_final_delta',...
       'MVMD_AC_reconstruction_relative','KF_seconds_PC',...
       'MVMD_seconds_PC','train_samples','validation_samples'});
    writetable(S,fullfile(out,'E05_solver_summary.csv'));
    % A/B/C use exactly the same fitted frequency and time indices.
    outputs=struct('raw',x,'KF_fixed',yKF);
    if filtReady,outputs.FVIM_MVMD_frozen=yMVMD;end
    save(fullfile(out,'E05_data.mat'),'outputs','va','tr','tv','cs','cfg','-v7');
    fig=figure('Visible','off');
    labels={'a / col4','s / col7'};
    for ch=1:2
        subplot(2,1,ch);
        plot(tv,x(:,ch),'Color',[.72 .72 .72]);hold on;
        plot(tv,yKF(:,ch),'b','LineWidth',1);
        if filtReady,plot(tv,yMVMD(:,ch),'r--','LineWidth',1);end
        grid on;ylabel(labels{ch});
        if ch==1,title([stem,' | same validation window']);end
        if ch==2,xlabel('Time (s)');end
        if filtReady,legend({'raw','KF fixed','MVMD frozen'},'Location','best');
        else,legend({'raw','KF fixed'},'Location','best');end
    end
    saveas(fig,fullfile(out,'E05_filtered_comparison.png'));close(fig);
    % Demodulation is only physically dimensional if each condition has
    % real, independently calibrated b and K values.
    calibrated=numel(cs.bias)==2 && numel(cs.scale)==2 && ...
      all(isfinite([cs.bias(:);cs.scale(:)])) && all(cs.scale>0);
    if calibrated
        phi=deg2rad(cs.phi0Deg)+cs.rollSign*2*pi*cs.fRoll*tv;
        dmNames=fieldnames(outputs);
        Tdem=table(tv,'VariableNames',{'time_s'});
        for j=1:numel(dmNames)
            name=dmNames{j}; q=outputs.(name);
            a=(q(:,1)-cs.bias(1))/cs.scale(1);
            s=(q(:,2)-cs.bias(2))/cs.scale(2);
            ox=s.*cos(phi)+a.*sin(phi);
            oy=s.*sin(phi)-a.*cos(phi);
            Tdem.([name '_omegaX_deg_s'])=ox;
            Tdem.([name '_omegaY_deg_s'])=oy;
            Tdem.([name '_omegaMagnitude_deg_s'])=hypot(ox,oy);
        end
        writetable(Tdem,fullfile(out,'E05_signed_transverse_rates.csv'));
        plotSigned(Tdem,out);
        % Only after explicit verified coordinate mapping.
        if cfg.reference.enabled
            E05_reference_attitude(M(va,1:3),tv,Tdem,cs,cfg.reference,out);
        else
            writeText(fullfile(out,'E05_attitude_status.txt'),...
              'SKIPPED orientation RMSE: gimbal axes, rotation signs, sensor installation, and Omega frame unverified. Signed OmegaX/Y still exported.');
        end
    else
        msg=sprintf(['UNCALIBRATED: %s\nSpecify per-condition measured bias=[ba bs] and scale=[Ka Ks] in E05_config.m.\n',...
            'Raw/KF/MVMD comparison is valid in sensor units, but deg/s and attitude metrics are not.\n'],stem);
        writeText(fullfile(out,'E05_REQUIRED_CALIBRATION.txt'),msg);
    end
    fprintf('Saved: %s\n',out);
end

function idx=clampIndices(N,fraction)
    idx=max(1,round(N*fraction(1))):min(N,round(N*fraction(2)));
    idx=idx(:);
end
function writeText(filename,msg)
    fid=fopen(filename,'w');if fid<0,error('Cannot write %s',filename);end
    fprintf(fid,'%s\n',msg);fclose(fid);
end
function plotSigned(T,out)
    names={'raw','KF_fixed','FVIM_MVMD_frozen'};
    fig=figure('Visible','off');
    for axis=1:2
       subplot(2,1,axis);hold on;
       for k=1:numel(names)
           c=[names{k} '_omegaX_deg_s'];
           if axis==2,c=[names{k} '_omegaY_deg_s'];end
           if any(strcmp(T.Properties.VariableNames,c))
               plot(T.time_s,T.(c));
           end
       end
       grid on;ylabel(sprintf('Omega_%s (deg/s)',char('x'+axis-1)));
       if axis==2,xlabel('Time(s)');end
       legend(names(1:(width(T)-1)/3),'Location','best');
    end
    saveas(fig,fullfile(out,'E05_signed_rates.png'));close(fig);
end
