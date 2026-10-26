%% E06-A | 双通道联合MVMD真实Z类动态边带保真检验
% 仅做信号层分析：不求Omega_x/Omega_y、不求姿态角、不重新标定、不做FVIM网格搜索。
% 基准：冻结E05 MVMD与CS/SampEn阈值；两个通道同时分解，模态分类成对执行。
clear; clc; close all;
cfg = E06_config();
if ~isfolder(cfg.dataRoot)
    error('E06:NoDataRoot', '请在 E06_config.m 修改 cfg.dataRoot: %s',cfg.dataRoot);
end
if ~isfolder(cfg.outputRoot), mkdir(cfg.outputRoot); end
allSummary = table();
for i = 1:numel(cfg.testFiles)
    filename = cfg.testFiles{i};
    fprintf('\n[E06-A %d/%d] %s\n',i,numel(cfg.testFiles),filename);
    try
        summary = runOneCase(filename,cfg);
        allSummary = [allSummary; summary]; %#ok<AGROW>
    catch ME
        fprintf(2,'E06-A FAILED %s: %s\n',filename,ME.message);
        [~,base] = fileparts(filename);
        caseFolder = fullfile(cfg.outputRoot,base);
        if ~isfolder(caseFolder), mkdir(caseFolder); end
        fid = fopen(fullfile(caseFolder,'ERROR.txt'),'w');
        if fid>0
            fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));
            fclose(fid);
        end
    end
end
if ~isempty(allSummary)
    writetable(allSummary,fullfile(cfg.outputRoot,'E06_run_summary.csv'));
end
fprintf('\nE06-A finished. Output: %s\n',cfg.outputRoot);

function summary = runOneCase(filename,cfg)
    tok = regexp(upper(filename),'^Z(\d+)-(\d+)\.TXT$','tokens','once');
    assert(~isempty(tok), 'E06-A只支持Zf-1.TXT、Zf-2.TXT类型');
    fRoll = str2double(tok{1});code=str2double(tok{2});
    assert(abs(fRoll-cfg.rollHz)<eps, '滚转频率不匹配：需要修改cfg.rollHz');
    % 下面的工况映射仅适用于本次目标文件，来自标定源码的工况说明。
    switch code
        case 1, fm=1; nominalAmplitudeDeg=1;
        case 2, fm=2; nominalAmplitudeDeg=3;
        otherwise, error('E06:UnknownMotionCode','此轮不分析Z%d-%d',fRoll,code);
    end
    src = findOneFile(cfg.dataRoot,filename);
    rec = read_gyro_txt(src);
    assert(rec.invalidRows==0,'存在被删除的无效行，可能影响时间轴');
    assert(rec.headerBlockRuleMatches,'TXT文件头块计数规则不满足');
    data = double(rec.data);
    first = 1 + round(cfg.skipSeconds*cfg.fs);
    last = min(size(data,1),first+round(cfg.maxDurationSeconds*cfg.fs)-1);
    assert(last-first+1 >= round(cfg.minDurationSeconds*cfg.fs), ...
        '有效记录不足%.1fs，请检查数据文件长度或调整skipSeconds',cfg.minDurationSeconds);
    Y = data(first:last,[4,7]);
    N=size(Y,1); t=((0:N-1)'/cfg.fs);
    outer = rad2deg(unwrap(deg2rad(data(first:last,1))));
    inner = rad2deg(unwrap(deg2rad(data(first:last,3))));
    [~,base] = fileparts(filename);
    folder = fullfile(cfg.outputRoot,base);
    if ~isfolder(folder),mkdir(folder);end
    targetF = [fRoll, abs(fRoll-fm), fRoll+fm];
    targetNames = {'roll_carrier','lower_sideband','upper_sideband'};
    assert(numel(unique(targetF))==numel(targetF), '目标频率重复，无法区分');
    % 外框角度独立检查，证明输入运动确实存在fm调制（仅作判据，不参加滤波）。
    outerH = harmonicLS(outer,t,[fm]);
    outerAmpDeg = abs(outerH.coeff(1));
    fprintf('  Outer-frame measured angle amplitude at %.3fHz: %.6f deg\n',fm,outerAmpDeg);
    tic;
    [modes,centers,info]=gyro_mvmd(Y,cfg.fs,cfg.mvmd.K,cfg.mvmd.alpha,...
        cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
    elapsed=toc;
    solver = table(string(base),N,fRoll,fm,nominalAmplitudeDeg,outerAmpDeg,...
        logical(info.converged),info.iterations,info.delta,info.AC_reconRelative,elapsed,...
        'VariableNames',{'caseName','samples','rollHz','motionHz',...
        'nominalOuterAngleAmplitudeDeg','measuredOuterAngleAmplitudeDeg',...
        'MVMD_converged','iterations','finalDelta','AC_reconstructionLoss','elapsed_s'});
    writetable(solver,fullfile(folder,'E06_solver.csv'));
    if ~info.converged
        fprintf(2,'  MVMD未收敛，结果不用于论文比较：iterations=%d\n',info.iterations);
        summary=solver;
        return;
    end
    rcfg=cfg;
    [Yf,details]=gyro_reconstruct(modes,centers,Y,rcfg,'paper');
    writetable(details,fullfile(folder,'E06_mode_details.csv'));
    modesTable=table();
    for k=1:cfg.mvmd.K
        for ch=1:2
            hm=harmonicLS(modes(:,ch,k),t,targetF);
            names={'a','s'};
            for j=1:numel(targetF)
                modesTable=[modesTable;table(k,string(names{ch}),centers(k),... 
                    details.class_0drop_1mix_2keep(k),string(targetNames{j}),targetF(j),...
                    abs(hm.coeff(j)),atan2d(imag(hm.coeff(j)),real(hm.coeff(j))),...
                    'VariableNames',{'IMF','channel','centerHz','class_0drop_1mix_2keep',...
                    'component','frequencyHz','harmonicAmplitude','coefficientPhaseDeg'})]; %#ok<AGROW>
            end
        end
    end
    writetable(modesTable,fullfile(folder,'E06_IMF_modes.csv'));
    rows=table();
    complexRaw=zeros(2,numel(targetF));complexFiltered=zeros(2,numel(targetF));
    signalStats=struct('raw',cell(1,2),'filtered',cell(1,2));
    channels={'a','s'};
    for ch=1:2
        r=harmonicLS(Y(:,ch),t,targetF);
        f=harmonicLS(Yf(:,ch),t,targetF);
        complexRaw(ch,:)=r.coeff;
        complexFiltered(ch,:)=f.coeff;
        signalStats(ch).raw=r;
        signalStats(ch).filtered=f;
        for j=1:numel(targetF)
            rawAmp=abs(r.coeff(j)); filtAmp=abs(f.coeff(j));
            retention=filtAmp/max(rawAmp,eps);
            phaseChange=wrap180(atan2d(imag(f.coeff(j)),real(f.coeff(j)))-...
                atan2d(imag(r.coeff(j)),real(r.coeff(j))));
            % 近似最小二乘系数SNR仅用于诊断，强相关噪声下会高估显著性。
            snrProxyDb=20*log10(max(rawAmp,eps)/max(r.stdCoeff,eps));
            passed=snrProxyDb>=cfg.sidebandReliabilityDb;
            attenuated=retention < (1-cfg.attenuationAlertPct/100);
            rows=[rows;table(string(base),string(channels{ch}),string(targetNames{j}),...
                targetF(j),rawAmp,filtAmp,100*retention,phaseChange,snrProxyDb,...
                passed,logical(passed && attenuated),...
                'VariableNames',{'caseName','channel','component','frequencyHz',...
                'rawAmplitude','filteredAmplitude','amplitudeRetentionPct',...
                'coefficientPhaseShiftDeg','rawHarmonicFitProxyDb',...
                'heuristicSidebandResolvable','attenuationOverThresholdIfResolvable'})]; %#ok<AGROW>
        end
    end
    writetable(rows,fullfile(folder,'E06_sideband_metrics.csv'));
    consistency=table();
    for j=1:numel(targetF)
        rawDiff=wrap180(atan2d(imag(complexRaw(2,j)),real(complexRaw(2,j)))-...
          atan2d(imag(complexRaw(1,j)),real(complexRaw(1,j))));
        filtDiff=wrap180(atan2d(imag(complexFiltered(2,j)),real(complexFiltered(2,j)))-...
          atan2d(imag(complexFiltered(1,j)),real(complexFiltered(1,j))));
        consistency=[consistency; table(string(base),string(targetNames{j}),targetF(j),...
            rawDiff,filtDiff,wrap180(filtDiff-rawDiff),...
            'VariableNames',{'caseName','component','frequencyHz','raw_S_minus_A_phaseDeg',...
            'filtered_S_minus_A_phaseDeg','relativePhaseChangeDeg'})]; %#ok<AGROW>
    end
    writetable(consistency,fullfile(folder,'E06_channel_consistency.csv'));
    noiseMetrics=table();
    for ch=1:2
        r=signalStats(ch).raw;f=signalStats(ch).filtered;
        % 非模型分量功率仅为噪声代理，动态有效信息也可能在拟合频率之外。
        rawResidual=sqrt(mean(r.residual.^2));
        filtResidual=sqrt(mean(f.residual.^2));
        diffRms=sqrt(mean((Yf(:,ch)-Y(:,ch)).^2));
        noiseMetrics=[noiseMetrics; table(string(base),string(channels{ch}),...
            rawResidual,filtResidual,20*log10(max(rawResidual,eps)/max(filtResidual,eps)),...
            diffRms,'VariableNames',{'caseName','channel',...
            'rawNonTargetResidualRMS_proxy','filteredNonTargetResidualRMS_proxy',...
            'nonTargetResidualSuppressionDb_proxy','rawMinusFilteredRMS_notTruthError'})]; %#ok<AGROW>
    end
    writetable(noiseMetrics,fullfile(folder,'E06_noise_metrics.csv'));
    % 图1：每路全部IMF在目标频段的幅值谱；仅展示，不参与分类。
    fig=figure('Visible','off','Color','w','Position',[100 100 1080 750]);
    for ch=1:2
        subplot(2,1,ch);
        [freq,rawSpect]=onesidedAmplitude(Y(:,ch),cfg.fs);
        plot(freq,rawSpect,'k','LineWidth',1.5,'DisplayName','raw');hold on;
        for k=1:cfg.mvmd.K
            [~,s]=onesidedAmplitude(modes(:,ch,k),cfg.fs);
            plot(freq,s,'LineWidth',0.8,'DisplayName',sprintf('IMF%d',k));
        end
        xlim([0 cfg.plotMaxHz]);grid on;
        ylabel(sprintf('Channel %s amplitude',channels{ch}));
        if ch==1,title(sprintf('%s | joint MVMD IMF spectra (0-%g Hz)',base,cfg.plotMaxHz));end
        for j=1:numel(targetF), xline(targetF(j),':','HandleVisibility','off');end
        if ch==1,legend('Location','eastoutside');end
    end
    xlabel('Frequency (Hz)');saveas(fig,fullfile(folder,'E06_IMF_spectrum.png'));close(fig);
    % 图2：原始与冻结MVMD输出幅值谱。
    fig=figure('Visible','off','Color','w','Position',[100 100 950 690]);
    for ch=1:2
        subplot(2,1,ch);
        [freq,s0]=onesidedAmplitude(Y(:,ch),cfg.fs);
        [~,s1]=onesidedAmplitude(Yf(:,ch),cfg.fs);
        plot(freq,s0,'Color',[.6 .6 .6],'DisplayName','raw');hold on;
        plot(freq,s1,'b','DisplayName','MVMD frozen');
        xlim([0 cfg.plotMaxHz]);grid on;
        ylabel(sprintf('%s amplitude',channels{ch}));
        for j=1:numel(targetF),xline(targetF(j),':','HandleVisibility','off');end
        if ch==1,title([base ' | raw vs MVMD']);legend('Location','northeast');end
    end
    xlabel('Frequency (Hz)');saveas(fig,fullfile(folder,'E06_spectrum_compare.png'));close(fig);
    % 图3：双通道部分时域数据。
    fig=figure('Visible','off','Color','w','Position',[100 100 930 670]);
    for ch=1:2
        subplot(2,1,ch);
        plot(t,Y(:,ch),'Color',[.65 .65 .65]);hold on;
        plot(t,Yf(:,ch),'b');grid on;
        ylabel(sprintf('%s raw units',channels{ch}));
        if ch==1,title([base ' | raw / MVMD frozen']);end
    end
    xlabel('time in selected block (s)');saveas(fig,fullfile(folder,'E06_waveform_compare.png'));close(fig);
    fid=fopen(fullfile(folder,'E06_interpretation.txt'),'w');
    if fid>0
        fprintf(fid,'REAL Z-class data: %s | roll=%gHz, motion=%gHz, outer encoder motion amplitude=%gdeg\n',filename,fRoll,fm,outerAmpDeg);
        fprintf(fid,'Frozen K=%d alpha=%.12g tau=%g; one JOINT MVMD decomposition; no FVIM optimizer.\n',cfg.mvmd.K,cfg.mvmd.alpha,cfg.mvmd.tau);
        fprintf(fid,'Expected candidate sidebands: %.6g and %.6gHz. Verify these with measured data; they are not guaranteed physical truth.\n',abs(fRoll-fm),fRoll+fm);
        fprintf(fid,'Raw-vs-filter harmonic retention uses simultaneous least-squares fitting of 3 frequencies. Raw is noisy reference, NOT noiseless truth.\n');
        fprintf(fid,'Coefficient phase convention: H=C+i*S in y=b+C*cos(2*pi*f*t)+S*sin(2*pi*f*t).\n');
        fprintf(fid,'Sideband resolvable flag is a heuristic single-tone fit proxy; autocorrelation may invalidate its nominal strength.\n');
        fprintf(fid,'The outer angle is only a PHYSICAL EXCITATION CHECK; it was not used for filtering, IMF selection or parameter fitting.\n');
        fprintf(fid,'Z-class dynamics are NOT already validated by static W-class calibration: no signed body-rate or attitude RMSE here.\n');
        fprintf(fid,'The same pre-frozen algorithm settings are used for Z3-1 and Z3-2. Do not tune on Z3-2 while calling it independent.\n');
        fclose(fid);
    end
    save(fullfile(folder,'E06_record.mat'),'Y','Yf','t','outer','inner','modes',...
        'centers','details','rows','consistency','noiseMetrics','info','cfg','-v7');
    fprintf('  converged at %d iterations; classes: %s\n',info.iterations,mat2str(details.class_0drop_1mix_2keep.'));
    fprintf('  Outputs: %s\n',folder);
    summary=solver;
end

function fit = harmonicLS(y,t,freqs)
    y=y(:);t=t(:);freqs=freqs(:).';
    assert(numel(y)==numel(t) && all(isfinite(y)),'非有限谐波拟合输入');
    X=ones(numel(t),1);
    for k=1:numel(freqs)
        X=[X cos(2*pi*freqs(k)*t) sin(2*pi*freqs(k)*t)]; %#ok<AGROW>
    end
    assert(rank(X)==size(X,2),'频率设计矩阵秩不足');
    beta=X\y;predicted=X*beta;
    coeff=complex(beta(2:2:end),beta(3:2:end)).';
    residual=y-predicted;
    % 白噪声OLS幅值的粗略系数标准差（仅代理）。
    sigma=sqrt(sum(residual.^2)/max(1,numel(y)-size(X,2)));
    sd=sigma * sqrt(2/numel(y));
    fit.coeff=coeff;fit.residual=residual;fit.stdCoeff=sd;
end

function [frequency,amplitude]=onesidedAmplitude(x,fs)
    x=x(:)-mean(x);
    N=numel(x);
    half=floor(N/2);
    F=fft(x);
    amplitude=(2/N)*abs(F(1:half+1));
    amplitude(1)=amplitude(1)/2;
    if mod(N,2)==0,amplitude(end)=amplitude(end)/2;end
    frequency=(0:half)'*fs/N;
end

function path=findOneFile(root,name)
    d=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'))];
    selected=arrayfun(@(x)strcmpi(x.name,name),d);
    d=d(selected);
    if isempty(d),error('E06:MissingFile','找不到%s，当前dataRoot=%s',name,root);end
    paths=fullfile({d.folder},{d.name});
    [~,ii]=unique(lower(string(paths)),'stable');d=d(ii);
    if numel(d)~=1,error('E06:DuplicateFile','发现多个同名%s，请处理',name);end
    path=fullfile(d.folder,d.name);
end
function a=wrap180(a)
    a=mod(a+180,360)-180;
end
