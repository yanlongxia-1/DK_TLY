%% E05 - calibration-aware REAL W data benchmark.
% Outputs: RAW/KF/MVMD carrier metrics, and signed EXTERNAL-GIMBAL
% effective angular-rate estimate. Not Omega_x/Omega_y body rates yet.
clear; clc; close all;
cfg=E05_calibrated_config();
if ~isfolder(cfg.dataRoot),error('Data folder not found: %s',cfg.dataRoot);end
if ~isfile(cfg.calibrationCsv),error('Calibration CSV not found: %s',cfg.calibrationCsv);end
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
cal=readtable(cfg.calibrationCsv,'TextType','string');
required={'operatingMode','rollHz','biasA','biasS','scaleACos','scaleASin',...
          'scaleSCos','scaleSSin','recommended','biasSourceFile'};
assert(all(ismember(required,cal.Properties.VariableNames)),...
    'The calibration table is missing required columns.');
for i=1:numel(cfg.testFiles)
    nm=char(cfg.testFiles{i});
    try
        oneCase(nm,cal,cfg);
    catch ME
        fprintf(2,'FAILED %s: %s\n',nm,ME.message);
        fid=fopen(fullfile(cfg.outputRoot,['ERROR_' regexprep(nm,'[^a-zA-Z0-9]','_') '.txt']),'w');
        if fid>0,fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));fclose(fid);end
    end
end
fprintf('All selected W cases processed. Results: %s\n',cfg.outputRoot);

function oneCase(name,cal,cfg)
    tok=regexp(upper(name),'^W(\d+)-(\d+)\.TXT$','tokens','once');
    if isempty(tok),error('Only Wf-rate.TXT constant-outer-rate files are supported now.');end
    fRoll=str2double(tok{1});
    mask=strcmpi(string(cal.operatingMode),'rolling') & (cal.rollHz==fRoll);
    assert(nnz(mask)==1,'Exactly one rolling calibration row required for %g Hz.',fRoll);
    c=cal(mask,:);
    assert(logical(c.recommended),'Calibration row not recommended; stop before quantitative reporting.');
    if ~all(isfinite([c.biasA,c.biasS,c.scaleACos,c.scaleASin,c.scaleSCos,c.scaleSSin]))
        error('Missing calibration coefficient.');
    end
    W=findOneFile(cfg.dataRoot,name);
    ename=char(c.biasSourceFile);
    E=findOneFile(cfg.dataRoot,ename);
    rec=read_gyro_txt(W); er=read_gyro_txt(E);
    assert(rec.invalidRows==0 && er.invalidRows==0,'Invalid rows may shift timestamps; stop.');
    if ~rec.headerBlockRuleMatches || ~er.headerBlockRuleMatches
        error('Unexpected header or valid-row count.');
    end
    y=double(rec.data(:,[4 7]));N=size(y,1);
    assert(N>=2000,'Too few data rows.');
    thetaW=rad2deg(unwrap(deg2rad(rec.data(:,3))));
    thetaE=er.data(:,3);
    % Fit the E-class carrier from INDEPENDENT ZERO-INPUT FILE, NOT from W test.
    [first,last]=trimIndices(size(er.data,1),cfg.fs,cfg);
    ha=robust_harmonic_fit(thetaE(first:last),er.data(first:last,4),cfg);
    hs=robust_harmonic_fit(thetaE(first:last),er.data(first:last,7),cfg);
    HE=[ha.cosCoefficient,ha.sinCoefficient;hs.cosCoefficient,hs.sinCoefficient];
    bE=[c.biasA c.biasS];
    EbiasFit=[ha.bias hs.bias];
    % Difference is a QC warning, not validation-time rebiasing.
    if max(abs(EbiasFit-bE))>0.005
        warning('%s: E-file fitted DC differs from CAL table by > 0.005 raw units.',name);
    end
    tr=fracIdx(N,cfg.trainFraction);
    va=fracIdx(N,cfg.validationFraction);
    assert(isempty(intersect(tr,va)),'Train and validation overlap.');
    raw=y(va,:); tAll=(0:N-1)'/cfg.fs;t=tAll(va);
    % Frozen oscillator KF: train-only R; MVMD once on validation for OFFLINE evaluation.
    tic; KF=E05_KF_fixed(y,cfg.fs,fRoll,tr,cfg.kf.qOscOverR,cfg.kf.qBiasOverR);
    elapsedKF=toc; KF=KF(va,:);
    tic;
    [modes,centers,info]=gyro_mvmd(raw,cfg.fs,cfg.mvmd.K,cfg.mvmd.alpha,...
       cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
    elapsedMVMD=toc;
    ok=logical(info.converged);
    [~,stem]=fileparts(name);
    out=fullfile(cfg.outputRoot,stem);
    if ~isfolder(out),mkdir(out);end
    results=struct('raw',raw,'KF_fixed',KF);
    if ok
        rcfg=cfg;rcfg.rollHz=fRoll;
        [v,details]=gyro_reconstruct(modes,centers,raw,rcfg,'paper');
        results.MVMD_frozen=v;
        writetable(details,fullfile(out,'E05_mode_details.csv'));
    end
    timeScore=table(string(stem),fRoll,N,numel(va),logical(info.converged),...
        info.iterations,info.delta,info.AC_reconRelative,elapsedKF,elapsedMVMD,...
        'VariableNames',{'case_name','roll_Hz','total_samples','evaluation_samples',...
        'MVMD_converged','MVMD_iterations','MVMD_delta','AC_recon_relative',...
        'KF_seconds','MVMD_seconds'});
    writetable(timeScore,fullfile(out,'E05_solver.csv'));
    methods=fieldnames(results);
    allMetrics=table();
    rateReport=table();
    thetaEncoder=thetaW(va);
    % Operational fixed-frequency phase. Assumes START angle and direction
    % are already provided by the bench for alignment; NO phase tracking.
    signRoll=sign(median(diff(thetaW))); assert(signRoll~=0,'Inner angle stationary');
    thetaFixed=thetaW(1)+signRoll*360*fRoll*t;
    % Angle encoder also gives independent outer GIMBAL speed ground truth.
    % For constant W inputs, estimate speed from robust slope across EVAL only
    % (EVALUATION LABEL, NEVER used to fit bias/scale/filter coefficients).
    outer=rad2deg(unwrap(deg2rad(rec.data(:,1))));
    outerFit=robust_irls([ones(numel(va),1) t],outer(va),cfg);
    omegaRef=outerFit.beta(2);
    for mi=1:numel(methods)
        method=methods{mi}; yy=results.(method);
        met=gyro_metrics(raw,yy,cfg.fs,fRoll,method);
        allMetrics=[allMetrics;met]; %#ok<AGROW>
        [omegaEnc,denEnc]=extractOuterRate(yy,thetaEncoder,HE,bE,c);
        [omegaFixed,denFixed]=extractOuterRate(yy,thetaFixed,HE,bE,c);
        % Reconstruction inputs for reference comparability are always the SAME.
        eEnc=sqrt(mean((omegaEnc-omegaRef).^2));
        eFixed=sqrt(mean((omegaFixed-omegaRef).^2));
        rateReport=[rateReport;table(string(method),omegaRef,mean(omegaEnc),...
           eEnc,mean(omegaFixed),eFixed,min(denEnc),min(denFixed),...
           'VariableNames',{'method','outer_encoder_rate_ref_dps',...
           'mean_estimate_encoder_phase_dps','RMSE_encoder_phase_deg_s',...
           'mean_estimate_fixed_phase_dps','RMSE_fixed_phase_deg_s',...
           'min_response_power_encoder','min_response_power_fixed'})]; %#ok<AGROW>
        if mi==1
            T=table(t,thetaEncoder,thetaFixed,repmat(omegaRef,numel(t),1),...
                'VariableNames',{'time_s','inner_angle_encoder_deg',...
                'inner_angle_fixed_frequency_deg','outer_encoder_rate_ref_deg_s'});
        end
        T.([method '_outer_equiv_encoder_deg_s'])=omegaEnc;
        T.([method '_outer_equiv_fixed_deg_s'])=omegaFixed;
    end
    writetable(allMetrics,fullfile(out,'E05_signal_metrics.csv'));
    writetable(rateReport,fullfile(out,'E05_outer_rate_metrics.csv'));
    writetable(T,fullfile(out,'E05_outer_rate_series.csv'));
    fig=figure('Visible','off');
    for j=1:2
        subplot(2,1,j);plot(t,raw(:,j),'Color',[.68 .68 .68]);hold on;
        plot(t,KF(:,j),'b');if ok,plot(t,results.MVMD_frozen(:,j),'r--');end
        if j==1,title(stem);end
        ylabel(sprintf('Channel %d raw units',j));grid on;
    end
    xlabel('t (s)');saveas(fig,fullfile(out,'E05_waveforms.png'));close(fig);
    fig=figure('Visible','off');hold on;
    for mi=1:numel(methods)
        method=methods{mi};plot(t,T.([method '_outer_equiv_encoder_deg_s']),'DisplayName',method);
    end
    yline(omegaRef,'k--','DisplayName','Outer gimbal encoder');
    legend('Location','best');xlabel('t (s)');ylabel('Omega_outer (deg/s)');grid on;
    title([stem,' | Encoder-referenced W-class effective input estimate']);
    saveas(fig,fullfile(out,'E05_outer_rate_compare.png'));close(fig);
    save(fullfile(out,'E05_record.mat'),'T','rateReport','allMetrics','results',...
        'c','HE','bE','va','tr','cfg','-v7');
    fid=fopen(fullfile(out,'E05_interpretation.txt'),'w');
    if fid>0
        fprintf(fid,['REAL W-class data. Frozen calibration trained with E + W10/30/90.\n',...
            'W20/W60 held out from CALIBRATION regression, NOT automatically independent of earlier FILTER exploration.\n',...
            'theta encoder phase: diagnostic upper-bound offline input; theta fixed: initialized from encoder start angle, constant roll frequency.\n',...
            'Omega_outer is an EQUIVALENT signed gimbal-axis input under the calibrated W model.\n',...
            'IT IS NOT general body Omega_x/Omega_y, and NOT true attitude RMSE.\n',...
            'Never use outer reference label to correct per-file biases or gains.\n',...
            'H=C+i*S is a COS/SIN coefficient pair: response(theta)=C*cos(theta)+S*sin(theta).\n',...
            'The calibration E-file nonzero harmonic H_E is subtracted before signed input inference.\n',...
            'MVMD is offline full-block; comparison is NOT proof of causal embedded readiness.\n']);
        fclose(fid);
    end
    fprintf('PASS %s: MVMD converged=%d | ref outer=%+.4f deg/s | results %s\n',...
        stem,ok,omegaRef,out);
end

function [gamma,den]=extractOuterRate(Y,thetaDeg,HE,bE,c)
    C=cosd(thetaDeg(:));S=sind(thetaDeg(:));
    zeroH=[C S]*HE.';
    % Order: a,s. Remove E-class DC and E-class periodic zero-input response.
    u=Y-repmat(bE,size(Y,1),1)-zeroH;
    Ga=C*c.scaleACos+S*c.scaleASin;
    Gs=C*c.scaleSCos+S*c.scaleSSin;
    den=Ga.^2+Gs.^2;
    if min(den)<1e-12,error('Degenerate calibrated W response.');end
    gamma=(Ga.*u(:,1)+Gs.*u(:,2))./den;
end

function path=findOneFile(root,name)
    dd=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'))];
    keep=arrayfun(@(a)strcmpi(a.name,name),dd);
    dd=dd(keep);
    if ~isempty(dd)
        paths=fullfile({dd.folder},{dd.name});
        lowered=cellfun(@lower,paths,'UniformOutput',false);
        [~,uniqueIndex]=unique(lowered,'stable');
        dd=dd(uniqueIndex);
    end
    if numel(dd)~=1,error('Expected exactly one file %s under %s, got %d',name,root,numel(dd));end
    path=fullfile(dd.folder,dd.name);
end
function idx=fracIdx(N,f)
    idx=(max(1,round(f(1)*N)):min(N,round(f(2)*N))).';
end
function [a,b]=trimIndices(N,fs,cfg)
    a=1+round(cfg.settlingSeconds*fs);b=N-round(cfg.endTrimSeconds*fs);
    if b-a+1<round(cfg.minimumUsableSeconds*fs),a=1;b=N;end
end
