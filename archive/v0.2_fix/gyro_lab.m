function gyro_lab(stage,cfg)
% One-line entry point used by E00...E06 scripts.
stage=upper(char(stage));
if ~exist(cfg.outputRoot,'dir'),mkdir(cfg.outputRoot);end
switch stage
    case 'E00', runE00(cfg);
    case 'E01', runE01(cfg);
    case 'E02', runMethod(cfg,'E02','paper');
    case 'E03', runMethod(cfg,'E03','physics');
    case 'E04', runE04(cfg);
    case 'E05', runE05(cfg);
    case 'E06', runE06(cfg);
    case 'E07', runE07(cfg);
    case 'E08', runE08(cfg);
    case 'E09', runE09(cfg);
    case 'E10', runE10(cfg);
    otherwise, error('Unsupported stage %s',stage);
end
end

function folder=resultDir(cfg, stage, filename)
[~,stem]=fileparts(filename);
folder=fullfile(cfg.outputRoot,stem,stage);
if ~exist(folder,'dir'),mkdir(folder);end
end

function runE00(cfg)
rng(cfg.seed,'twister');
root=fullfile(cfg.root,'data');
if ~exist(root,'dir'),mkdir(root);end
freqs=[3,5,10,19];
for ff=freqs
    f=ff; fs=cfg.fs; N=8000;
    t=(0:N-1)'/fs;
    phi=2*pi*f*t+0.08*sin(2*pi*0.31*t);
    ox=20+5*sin(2*pi*0.54*t);
    oy=8+3*cos(2*pi*0.67*t);
    K=[0.0045,0.0050]; bias=[-0.186041,-0.072559];
    s=K(2)*(ox.*cos(phi)+oy.*sin(phi));
    a=K(1)*(ox.*sin(phi)-oy.*cos(phi));
    clean=[a,s];
    noisy=clean+bias+0.014*randn(N,2);
    spikeIdx=round(linspace(500,N-500,6));
    noisy(spikeIdx,1)=noisy(spikeIdx,1)+0.04;
    noisy(spikeIdx+11,2)=noisy(spikeIdx+11,2)-0.04;
    A=zeros(N,7);
    A(:,4)=noisy(:,1); A(:,7)=noisy(:,2);
    if f==3
        stem='synthetic_dualgyro';
    else
        stem=sprintf('synthetic_dualgyro_%dHz',f);
    end
    fpath=fullfile(root,[stem,'.txt']);
    writematrix(A,fpath,'Delimiter','tab');
    truth.fs=fs;truth.f0=f;truth.phi=phi;truth.omegaX=ox;truth.omegaY=oy;
    truth.clean=clean;truth.K=K;truth.bias=bias;
    save(fullfile(root,[stem,'_truth.mat']),'truth');
    if f==3
        sNorm=s/K(2); aNorm=a/K(1);
        oxRec=sNorm.*cos(phi)+aNorm.*sin(phi);
        oyRec=sNorm.*sin(phi)-aNorm.*cos(phi);
        maxerr=max(abs([oxRec-ox;oyRec-oy]),[],'all');
        fprintf('E00 demod exact-physics check max abs err: %.6g deg/s\n',maxerr);
        if maxerr>1e-9,error('Synthetic demod inversion failed');end
        out=resultDir(cfg,'E00',fpath);
        fid=fopen(fullfile(out,'selftest.txt'),'w');
        fprintf(fid,'Analytic inverse check maximum absolute error = %.14g deg/s\n',maxerr);
        fprintf(fid,'Truth phi includes small phase modulation; fixed f*t is deliberately imperfect.\n');
        fclose(fid);
        % Verify that the numerical MVMD/IMF pipeline runs end-to-end.
        smoke=noisy(1:2048,:);
        [sm,sc,si]=gyro_mvmd(smoke,fs,3,1000,cfg.mvmd.finalIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
        if ~si.converged
            error('E00: MVMD selftest failed convergence: delta=%g after %d iterations',si.delta,si.iterations);
        end
        localCfg=cfg;localCfg.rollHz=3;
        [sr,st]=gyro_reconstruct(sm,sc,smoke,localCfg,'paper');
        if any(~isfinite(sr(:))) || ~isequal(size(sr),size(smoke))
            error('MVMD smoke test failed: non-finite or wrong-size result');
        end
        writetable(st,fullfile(out,'selftest_imf.csv'));
        writetable(gyro_metrics(smoke,sr,fs,3,'smoke_test'),...
            fullfile(out,'selftest_metrics.csv'));
        fid=fopen(fullfile(out,'selftest.txt'),'a');
        fprintf(fid,'MVMD smoke test: iterations=%d; converged=%d; AC reconstruction residual=%g; final delta=%g\n',...
            si.iterations,si.converged,si.AC_reconRelative,si.delta);
        fclose(fid);
        fig=figure('Visible','off');
        plot(t(1:2000),noisy(1:2000,1),'-',t(1:2000),noisy(1:2000,2),'-');
        xlabel('Time (s)');ylabel('Sensor output (arb.)');legend('a / col4','s / col7');
        title('Synthetic dual-orthogonal gyro, raw including bias');grid on;
        saveas(fig,fullfile(out,'synthetic_preview.png'));close(fig);
    end
    fprintf('Created %s\n',fpath);
end
% 3Hz alias used in the multi-condition test.
src=fullfile(root,'synthetic_dualgyro.txt');
dst=fullfile(root,'synthetic_dualgyro_3Hz.txt');
copyfile(src,dst);
copyfile(fullfile(root,'synthetic_dualgyro_truth.mat'),fullfile(root,'synthetic_dualgyro_3Hz_truth.mat'));
fprintf('E00 completed. Run E01, then E02.\n');
end

function runE01(cfg)
d=gyro_read(cfg.file,cfg);out=resultDir(cfg,'E01',cfg.file);
raw=gyro_sinefit(d.y,cfg.fs,cfg.rollHz);
f0estimate=dominantFreq(d.y,cfg.fs);
ph=rad2deg(angle(exp(1i*(raw.phase(2)-raw.phase(1)))));
C=table({'a_col4';'s_col7'},mean(d.y,1)',std(d.y,0,1)',... 
    median(d.y,1)',min(d.y,[],1)',max(d.y,[],1)',...
    raw.amplitude(:),rad2deg(raw.phase(:)),raw.snrProxy(:),...
    'VariableNames',{'channel','mean_not_calibrated_bias','std','median','min','max',...
    'fit_amp_at_cfg_rollHz','fit_phase_deg','SNR_proxy_dB'});
writetable(C,fullfile(out,'data_audit.csv'));
info=table(d.N, d.N/cfg.fs, d.removedRows, cfg.fs, cfg.rollHz, f0estimate,ph,...
    'VariableNames',{'nSamples','duration_s','removedRows','fs_assumed','f0_cfg','joint_peak_freq_est','phase_s_minus_a_deg'});
writetable(info,fullfile(out,'data_audit_summary.csv'));
% Extra observable diagnostics; these cannot prove timestamp continuity.
dy=diff(d.y,1,1);sigD=1.4826*median(abs(dy-median(dy,1)),1);
spikes=sum(abs(dy)>8*(sigD+eps),1);
repeats=sum(abs(dy)<1e-12,1);
quality=table({'a_col4';'s_col7'},spikes(:),repeats(:),...
    'VariableNames',{'channel','large_first_difference_count','identical_adjacent_count'});
writetable(quality,fullfile(out,'data_quality_flags.csv'));
if abs(f0estimate-cfg.rollHz)>0.15
    warning('Estimated peak %.3fHz differs from configured roll %.3fHz; inspect before E02!',...
        f0estimate,cfg.rollHz);
end
fig=figure('Visible','off');
subplot(2,1,1);plot(d.t,d.y);grid on;legend('a','s');ylabel('output');title('Raw two-channel time series');
subplot(2,1,2);N=d.N;f=(0:floor(N/2))'*cfg.fs/N;
V=abs(fft(d.y-mean(d.y,1),[],1));plot(f,V(1:numel(f),:));
xlim([0,min(60,cfg.fs/2)]);grid on;ylabel('FFT magnitude');xlabel('Frequency (Hz)');
saveas(fig,fullfile(out,'raw_wave_spectrum.png'));close(fig);
fprintf('E01: N=%d, duration=%.3fs, peak freq ~ %.5fHz, phase s-a at cfg.rollHz = %.3f deg\n',...
    d.N,d.N/cfg.fs,f0estimate,ph);
fprintf('NOTICE: CSV mean is NOT a verified bias; sample continuity cannot be inferred from historical file without actual timestamps.\n');
fprintf('Outputs: %s\n',out);
end

function f0=dominantFreq(y,fs)
N=size(y,1);w=0.5-0.5*cos(2*pi*(0:N-1)'/(N-1));
F=sum(abs(fft((y-mean(y,1)).*w,[],1)).^2,2);
lo=max(2,ceil(0.5*N/fs)+1);hi=min(floor(N/2),floor(60*N/fs)+1);
[~,p]=max(F(lo:hi));idx=lo+p-1;
if idx>1 && idx<length(F)
    L=log(F(idx-1:idx+1)+eps);
    delta=0.5*(L(1)-L(3))/(L(1)-2*L(2)+L(3)+eps);
    delta=max(-0.5,min(0.5,delta));
else,delta=0;end
f0=(idx-1+delta)*fs/N;
end

function runMethod(cfg,stage,strategy)
d=gyro_read(cfg.file,cfg);
[trainIdx,valIdx]=splitIndices(d.N,cfg);
out=resultDir(cfg,stage,cfg.file);
xTrain=d.y(trainIdx,:); xVal=d.y(valIdx,:);
fprintf('%s: training %d points, validation %d points, optimizer %s\n',...
    stage,numel(trainIdx),numel(valIdx),strategy);
timer=tic;
[k,alpha,trace,fitinfo]=gyro_fvim(xTrain,cfg,strategy);
searchSec=toc(timer);
timer=tic;
[modes,centers,mvinfo]=gyro_mvmd(xVal,cfg.fs,k,alpha,...
    cfg.mvmd.finalIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
if ~mvinfo.converged
    error('E02/E03:nonconvergence', ...
        'MVMD validation failed convergence: delta=%g, centerDeltaHz=%g, iter=%d. Results NOT saved.', ...
        mvinfo.delta,mvinfo.centerDeltaHz,mvinfo.iterations);
end
[yFilt,imfTable]=gyro_reconstruct(modes,centers,xVal,cfg,strategy);
filterSec=toc(timer);
metrics=gyro_metrics(xVal,yFilt,cfg.fs,cfg.rollHz,stage);
writetable(metrics,fullfile(out,'filter_metrics.csv'));
writetable(imfTable,fullfile(out,'imf_classification.csv'));
searchHistory=array2table(trace,'VariableNames',{'iteration','K','alpha','fitness'});
writetable(searchHistory,fullfile(out,'fvim_trace.csv'));
summary=table(k,alpha,searchSec,filterSec,fitinfo.calls,mvinfo.iterations, ...
    mvinfo.converged,mvinfo.delta,mvinfo.centerDeltaHz,mvinfo.AC_reconRelative, ...
    mvinfo.spectralResidual,fitinfo.boundaryK, ...
    sum(imfTable.class_0drop_1mix_2keep==0), ...
    sum(imfTable.class_0drop_1mix_2keep==1), ...
    sum(imfTable.class_0drop_1mix_2keep==2), ...
    'VariableNames',{'K_opt','alpha_opt','search_seconds_PC','filter_seconds_PC',...
    'fitness_calls','mvmd_iterations','mvmd_converged','mvmd_final_delta', ...
    'mvmd_center_change_Hz','MVMD_AC_reconstruction_relative', ...
    'MVMD_spectral_reconstruction_relative','FVIM_K_on_search_boundary', ...
    'n_discarded','n_mixed','n_kept'});
writetable(summary,fullfile(out,'method_summary.csv'));
% Indices are persisted to prevent accidental train/validation confusion.
record.filename=cfg.file; record.valIdx=valIdx;record.trainIdx=trainIdx;
record.raw=xVal;record.filt=yFilt;record.modes=modes;record.centers=centers;
record.k=k;record.alpha=alpha;record.strategy=strategy;record.cfg=cfg;
record.mvmdInfo=mvinfo; record.version='v0.2';
save(fullfile(out,'record.mat'),'record','metrics','summary','imfTable','-v7.3');
writematrix([d.t(valIdx),xVal,yFilt],fullfile(out,'validation_raw_filtered.csv'));
plotMethod(d.t(valIdx),xVal,yFilt,modes,centers,stage,out);
fprintf('%s: K=%d, alpha=%.3f, MVMD CONVERGED=%d (iter=%d, ACerr=%.3f%%)\n',...
    stage,k,alpha,mvinfo.converged,mvinfo.iterations,100*mvinfo.AC_reconRelative);
if fitinfo.boundaryK
    warning('FVIM optimum K lies on allowed boundary; E02B grid and multi-seed checks required before final claims.');
end
disp(metrics(:,{'channel','carrier_amp_keep_pct','phase_shift_deg','snr_proxy_gain_dB','relative_phase_change_deg'}));
fprintf('Outputs: %s\n',out);
end

function [tr,va]=splitIndices(N,cfg)
tr=max(1,round(N*cfg.trainFraction(1))):min(N,round(N*cfg.trainFraction(2)));
va=max(1,round(N*cfg.valFraction(1))):min(N,round(N*cfg.valFraction(2)));
if ~isempty(intersect(tr,va)),error('Train/test data overlap');end
if numel(tr)<800 || numel(va)<800
    warning('Less than 800 train or test samples; increase recording duration.');
end
end

function plotMethod(t,r,y,modes,centers,stage,out)
fig=figure('Visible','off');
subplot(2,1,1);plot(t,r(:,1),'-',t,y(:,1),'--');grid on;
title([stage,' | a / col4']);legend('measured','filtered');ylabel('signal');
subplot(2,1,2);plot(t,r(:,2),'-',t,y(:,2),'--');grid on;
title([stage,' | s / col7']);legend('measured','filtered');xlabel('Time (s)');ylabel('signal');
saveas(fig,fullfile(out,'raw_filtered.png'));close(fig);
K=size(modes,3);fig=figure('Visible','off');
for k=1:K
    subplot(K,1,k);plot(t,modes(:,1,k));grid on;
    ylabel(sprintf('IMF%d',k));
    title(sprintf('Joint MVMD IMF %d center %.2f Hz (a-channel shown)',k,centers(k)));
end
xlabel('Time(s)');saveas(fig,fullfile(out,'mvmd_imfs.png'));close(fig);
end

function runE04(cfg)
files=cfg.multiFiles;
for i=1:numel(files)
    if ~isfile(files{i})
        fprintf('E04 missing: %s\nRun E00 or edit cfg.multiFiles.\n',files{i});return;
    end
end
if isempty(files),error('cfg.multiFiles is empty');end
originCfg=cfg;originCfg.file=files{1};
first=gyro_read(files{1},originCfg);
[tr,~]=splitIndices(first.N,cfg);
f0=frequencyFromName(files{1},cfg.rollHz);
originCfg.rollHz=f0;
fprintf('E04 fitting both methods ONLY on first file %s\n',files{1});
[kb,ab]=gyro_fvim(first.y(tr,:),originCfg,'paper');
[kp,ap]=gyro_fvim(first.y(tr,:),originCfg,'physics');
rows={};
for i=1:numel(files)
    c=cfg;c.file=files{i};c.rollHz=frequencyFromName(files{i},cfg.rollHz);
    d=gyro_read(c.file,c);[~,idx]=splitIndices(d.N,c);xx=d.y(idx,:);
    for mm=1:2
        if mm==1,mode='paper';k=kb;alpha=ab; else,mode='physics';k=kp;alpha=ap;end
        [m,cent,~]=gyro_mvmd(xx,c.fs,k,alpha,c.mvmd.finalIterations,c.mvmd.tolerance,c.mvmd.tau);
        yy=gyro_reconstruct(m,cent,xx,c,mode);
        T=gyro_metrics(xx,yy,c.fs,c.rollHz,mode);
        for j=1:2
            rows(end+1,:)={files{i},c.rollHz,mode,T.channel{j},k,alpha,...
                T.carrier_amp_keep_pct(j),T.phase_shift_deg(j),...
                T.relative_phase_change_deg(j),T.snr_proxy_gain_dB(j)}; %#ok<AGROW>
        end
    end
end
out=resultDir(cfg,'E04',files{1});
Tab=cell2table(rows,'VariableNames',{'file','rollHz','method','channel',...
    'K_fixed','alpha_fixed','amplitude_keep_pct','phase_deg','relative_phase_change_deg','snr_proxy_gain_dB'});
writetable(Tab,fullfile(out,'cross_condition_fixed_parameter.csv'));
fprintf('E04 completed. Parameters trained only on first file; not reoptimized on test files.\n');
fprintf('Outputs: %s\n',out);
end

function f=frequencyFromName(filename,def)
[~,stem]=fileparts(filename);
tokens=regexp(stem,'_(\d+(?:\.\d+)?)Hz','tokens','once');
if isempty(tokens),f=def;else,f=str2double(tokens{1});end
end

function runE05(cfg)
base02=fullfile(resultDir(cfg,'E02',cfg.file),'record.mat');
base03=fullfile(resultDir(cfg,'E03',cfg.file),'record.mat');
if ~isfile(base02)||~isfile(base03)
    error('Run E02 and E03 for this SAME file before E05.');
end
b=load(base02,'record');p=load(base03,'record');
r=b.record;q=p.record;
if numel(r.valIdx)~=numel(q.valIdx)||any(r.valIdx~=q.valIdx)
    error('E02 and E03 validation indices mismatch');
end
[folder,stem]=fileparts(cfg.file);
truthFile=fullfile(folder,[stem,'_truth.mat']);
out=resultDir(cfg,'E05',cfg.file);
if isfile(truthFile)
    temp=load(truthFile,'truth');truth=temp.truth;
    k=truth.K; bias=truth.bias;
    if strcmpi(cfg.biasMode,'manual')
        if any(abs(cfg.bias(:)-bias(:))>1e-7)
            warning('cfg.bias differs from simulator truth.bias; check calibration');
        end
        subtract=[0,0];
    else
        subtract=bias;
    end
    phis={truth.phi(r.valIdx),'reference_phi'; ...
        2*pi*cfg.rollHz*((r.valIdx(:)-1)/cfg.fs)+cfg.phi0,'fixed_frequency_phi'};
    referencePresent=true;
elseif all(isfinite(cfg.sensitivity)) && strcmpi(cfg.biasMode,'manual')
    k=cfg.sensitivity;subtract=[0,0];
    phis={2*pi*cfg.rollHz*((r.valIdx(:)-1)/cfg.fs)+cfg.phi0,'assumed_phi'};
    referencePresent=false;
else
    fprintf('E05: Missing calibrated Ka/Ks AND verified offset and roll phase.\n');
    fprintf('No deg/s or angle RMSE can be asserted. Please configure sensitivity and bias, then rerun.\n');
    fid=fopen(fullfile(out,'REQUIRED_CALIBRATION.txt'),'w');
    fprintf(fid,'Required: verified Ka,Ks; zero offsets; roll phase vs time; coordinate/frame truth for angle RMSE.\n');
    fclose(fid);return;
end
methods={'raw','E02_paper','E03_physics'};
yset={r.raw,r.filt,q.filt}; rows={};
fig=figure('Visible','off');
for pp=1:size(phis,1)
    phi=phis{pp,1};label=phis{pp,2};
    for j=1:3
        v=yset{j}-subtract;
        a=v(:,1)/k(1);s=v(:,2)/k(2);
        ox=s.*cos(phi)+a.*sin(phi);
        oy=s.*sin(phi)-a.*cos(phi);
        mag=hypot(ox,oy);
        if referencePresent
            dx=ox-truth.omegaX(r.valIdx);
            dy=oy-truth.omegaY(r.valIdx);
            dm=mag-hypot(truth.omegaX(r.valIdx),truth.omegaY(r.valIdx));
            rx=sqrt(mean(dx.^2)); ry=sqrt(mean(dy.^2));rm=sqrt(mean(dm.^2));
        else
            rx=NaN;ry=NaN;rm=NaN;
        end
        rows(end+1,:)={methods{j},label,rx,ry,rm,...
            mean(ox),mean(oy),mean(mag),numel(ox)}; %#ok<AGROW>
        if pp==1
            subplot(2,1,1);hold on;plot(ox);ylabel('Omega_x (deg/s)');grid on;
            subplot(2,1,2);hold on;plot(oy);ylabel('Omega_y (deg/s)');grid on;
        end
    end
end
subplot(2,1,1);legend(methods,'Location','best');
subplot(2,1,2);legend(methods,'Location','best');xlabel('validation sample');
saveas(fig,fullfile(out,'demodulated_transverse_rates.png'));close(fig);
R=cell2table(rows,'VariableNames',{'method','phi_source',...
    'OmegaX_RMSE_deg_s','OmegaY_RMSE_deg_s','OmegaMagnitude_RMSE_deg_s',...
    'OmegaX_mean','OmegaY_mean','magnitude_mean','nSamples'});
writetable(R,fullfile(out,'demodulation_results.csv'));
if ~referencePresent
    warning('Output rates use user-calibrated constants and assumed phi; NO real ground-truth RMSE.');
end
fprintf('E05 completed. These are TRANSVERSE-RATE results, NOT Euler-angle RMSE.\n');
disp(R);
fprintf('Outputs: %s\n',out);
end

function runE06(cfg)
path=fullfile(resultDir(cfg,'E02',cfg.file),'record.mat');
if ~isfile(path),error('Run E02 first');end
S=load(path,'record'); R=S.record;
d=gyro_read(cfg.file,cfg);
windows=[1024,2048,4096];
rows=[];
for N=windows
    if d.N<N,continue;end
    x=d.y(1:N,:);
    for rep=1:2
        tic;
        [m,cent,infostr]=gyro_mvmd(x,cfg.fs,R.k,R.alpha,cfg.mvmd.finalIterations,...
            cfg.mvmd.tolerance,cfg.mvmd.tau);
        gyro_reconstruct(m,cent,x,cfg,'paper');
        sec=toc;
        rows(end+1,:)=[N,rep,sec,N/cfg.fs,sec/(N/cfg.fs),...
            infostr.iterations,double(infostr.converged)]; %#ok<AGROW>
    end
end
out=resultDir(cfg,'E06',cfg.file);
T=array2table(rows,'VariableNames',{'window_samples','repeat','PC_processing_seconds',...
    'data_window_seconds','compute_over_window_ratio','mvmd_iterations','converged'});
writetable(T,fullfile(out,'window_runtime.csv'));
fprintf('E06: offline PC timing; this does NOT establish an STM32H723 throughput bound.\n');
disp(T);
fprintf('Outputs: %s\n',out);
end


function runE07(cfg)
% 2x2 ablation: parameter objective versus mode-classification protection.
p2=fullfile(resultDir(cfg,'E02',cfg.file),'record.mat');
p3=fullfile(resultDir(cfg,'E03',cfg.file),'record.mat');
if ~isfile(p2)||~isfile(p3), error('Run E02 and E03 before E07');end
s2=load(p2,'record');s3=load(p3,'record');
r=s2.record;z=s3.record; rows={};
for i=1:2
    if i==1,k=r.k;alpha=r.alpha;objective='paper';
    else,k=z.k;alpha=z.alpha;objective='physics';end
    [m,cent,~]=gyro_mvmd(r.raw,cfg.fs,k,alpha,cfg.mvmd.finalIterations,...
        cfg.mvmd.tolerance,cfg.mvmd.tau);
    for j=1:2
        if j==1,selector='paper';else,selector='physics';end
        [yy,~]=gyro_reconstruct(m,cent,r.raw,cfg,selector);
        tag=[objective,'_opt__',selector,'_select'];
        T=gyro_metrics(r.raw,yy,cfg.fs,cfg.rollHz,tag);
        for cc=1:2
            rows(end+1,:)={objective,selector,T.channel{cc},k,alpha,...
                T.carrier_amp_keep_pct(cc),T.phase_shift_deg(cc),...
                T.relative_phase_change_deg(cc),T.snr_proxy_gain_dB(cc),...
                T.RMSE_vs_measured(cc)}; %#ok<AGROW>
        end
    end
end
out=resultDir(cfg,'E07',cfg.file);
Tab=cell2table(rows,'VariableNames',{'objective','IMF_classification','channel','K',...
    'alpha','amplitude_keep_pct','phase_deg','relative_phase_change_deg',...
    'snr_proxy_gain_dB','RMSE_vs_measured'});
writetable(Tab,fullfile(out,'two_by_two_ablation.csv'));
disp(Tab);
fprintf('E07: separates physical FVIM optimization and physical mode protection.\n');
fprintf('Outputs: %s\n',out);
end

function runE08(cfg)
% TRUE ground truth simulated noise stress. Fixed learned parameters only.
[folder,stem]=fileparts(cfg.file);
truthFile=fullfile(folder,[stem,'_truth.mat']);
if ~isfile(truthFile)
    fprintf('E08 requires synthetic truth from E00. Run with a synthetic file.\n');return;
end
p2=fullfile(resultDir(cfg,'E02',cfg.file),'record.mat');
p3=fullfile(resultDir(cfg,'E03',cfg.file),'record.mat');
if ~isfile(p2)||~isfile(p3),error('Run E02 and E03 first');end
S=load(truthFile,'truth');truth=S.truth;
S2=load(p2,'record');S3=load(p3,'record');
B=S2.record;P=S3.record;
[~,idx]=splitIndices(size(truth.clean,1),cfg);
idx=idx(1:min(2048,numel(idx)));
clean=truth.clean(idx,:);
noiseSigma=[0.005,0.014,0.035];
rows={};rng(cfg.seed+88,'twister');
for j=1:numel(noiseSigma)
    for rep=1:3
        sigma=noiseSigma(j);
        y=clean+truth.bias+sigma*randn(size(clean));
        if rep==3
            pos=round(linspace(80,size(y,1)-80,5));
            y(pos,1)=y(pos,1)+6*sigma;
            y(pos+5,2)=y(pos+5,2)-6*sigma;
        end
        for method=1:3
            if method==1
                tag='raw';v=y;
            elseif method==2
                tag='paper';[m,c,~]=gyro_mvmd(y,cfg.fs,B.k,B.alpha,...
                    cfg.mvmd.finalIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
                v=gyro_reconstruct(m,c,y,cfg,'paper');
            else
                tag='physics';[m,c,~]=gyro_mvmd(y,cfg.fs,P.k,P.alpha,...
                    cfg.mvmd.finalIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
                v=gyro_reconstruct(m,c,y,cfg,'physics');
            end
            % Compare true CLEAN SENSOR SIGNAL after removing known bias.
            err=v-truth.bias-clean;
            rmseA=sqrt(mean(err(:,1).^2));
            rmseS=sqrt(mean(err(:,2).^2));
            fitC=gyro_sinefit(clean,cfg.fs,cfg.rollHz);
            fitV=gyro_sinefit(v-truth.bias,cfg.fs,cfg.rollHz);
            amp=max(abs(fitV.amplitude./(fitC.amplitude+eps)-1));
            delta=rad2deg(angle(exp(1i*((fitV.phase(2)-fitV.phase(1))...
                -(fitC.phase(2)-fitC.phase(1))))));
            rows(end+1,:)={sigma,rep,tag,rmseA,rmseS,amp*100,delta}; %#ok<AGROW>
        end
    end
end
out=resultDir(cfg,'E08',cfg.file);
T=cell2table(rows,'VariableNames',{'noise_sigma','repeat','method',...
    'TRUE_sensor_RMSE_a','TRUE_sensor_RMSE_s',...
    'TRUE_max_carrier_amp_error_pct','TRUE_relative_phase_error_deg'});
writetable(T,fullfile(out,'synthetic_noise_robustness.csv'));
disp(T);
fprintf('E08: TRUE sensor reference used; not a real turntable angle accuracy claim.\n');
fprintf('Outputs: %s\n',out);
end


function runE09(cfg)
% Re-implemented oscillator KF comparator for independent repeatability.
% The user's original tuned KF must be separately imported to claim exact parity.
d=gyro_read(cfg.file,cfg);
[tr,va]=splitIndices(d.N,cfg);
x=d.y(va,:);
R=zeros(1,2);
for cc=1:2
    dx=diff(d.y(tr,cc));
    R(cc)=max(median(abs(dx-median(dx)))^2/(0.67448975^2*2),eps);
end
qOscRatio=1e-6; qBiasRatio=1e-8;
timer=tic;
y=gyro_kf(x,[cfg.fs,cfg.rollHz],R,qOscRatio,qBiasRatio);
sec=toc(timer);
T=gyro_metrics(x,y,cfg.fs,cfg.rollHz,'KF_independent_reimplementation');
out=resultDir(cfg,'E09',cfg.file);
writetable(T,fullfile(out,'KF_metrics.csv'));
writematrix([d.t(va),x,y],fullfile(out,'KF_raw_filtered.csv'));
parameters=table(R(1),R(2),qOscRatio,qBiasRatio,sec, ...
    'VariableNames',{'R_a_est','R_s_est','qOsc_over_R','qBias_over_R','PC_time_s'});
writetable(parameters,fullfile(out,'KF_parameters.csv'));
fprintf('E09: independent model KF, not exact reconstruction of old tuned code.\n');
disp(T);
fprintf('Outputs: %s\n',out);
end


function runE10(cfg)
% Allan deviation of stationary raw measurements; no false dynamic-run ARW.
if ~cfg.staticData
    error('E10 requires verified STATIONARY recording. Set cfg.staticData=true only after verification.');
end
d=gyro_read(cfg.file,cfg);
if d.N/cfg.fs < cfg.minAllanSeconds
    error('E10 needs >=%.0f sec stationary series, current duration %.2f sec.',...
        cfg.minAllanSeconds,d.N/cfg.fs);
end
if all(isfinite(cfg.sensitivity)) && all(cfg.sensitivity~=0)
    y=d.y./reshape(cfg.sensitivity,1,2);
    unit='deg_per_s';
else
    y=d.y;
    unit='sensor_units';
    warning('Sensitivity uncalibrated. Allan deviation is NOT in deg/s or deg/sqrt(h).');
end
N=size(y,1);
maxM=floor(N/20);
ms=unique(max(1,round(logspace(0,log10(maxM),26))));
allan=zeros(numel(ms),2);
for jj=1:numel(ms)
    m=ms(jj);
    nGroups=floor(N/m);
    yy=reshape(y(1:nGroups*m,:),m,nGroups,2);
    avg=reshape(mean(yy,1),nGroups,2);
    allan(jj,:)=sqrt(0.5*mean(diff(avg,1,1).^2,1));
end
tau=ms(:)/cfg.fs;
out=resultDir(cfg,'E10',cfg.file);
T=table(tau,allan(:,1),allan(:,2),...
    'VariableNames',{'tau_sec','AllanDeviation_a','AllanDeviation_s'});
writetable(T,fullfile(out,'static_Allan_deviation.csv'));
fig=figure('Visible','off');
loglog(tau,allan(:,1),'-o',tau,allan(:,2),'-s');grid on;
xlabel('Cluster averaging time tau (s)');
ylabel(['Allan deviation (',strrep(unit,'_',' '),')']);
legend('a','s','Location','best');
title('Static, nonoverlapping Allan deviation (NOT automatic ARW fit)');
saveas(fig,fullfile(out,'static_Allan_deviation.png'));close(fig);
fid=fopen(fullfile(out,'unit_and_limits.txt'),'w');
fprintf(fid,'unit=%s\n',unit);
fprintf(fid,'staticData flag=true; operator must verify actual stationary condition.\n');
fprintf(fid,'Non-overlapping Allan estimator; this file contains no automatic ARW coefficient.\n');
fprintf(fid,'Do not interpret brief, dynamic roll data as stationary random-walk performance.\n');
fclose(fid);
fprintf('E10 outputs: %s\n',out);
end
