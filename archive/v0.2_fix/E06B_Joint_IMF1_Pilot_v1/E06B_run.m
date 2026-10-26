function E06B_run()
%E06B_RUN 用E06-A保存的真实联合MVMD模态直接进行单次改进试验。
% E06B不再执行MVMD、不再跑FVIM；严格先冻结方案再对两组Z数据作结果评价。
clc;cfg=E06B_config();
if ~isfolder(cfg.outputRoot),mkdir(cfg.outputRoot);end
E06B_selftest();
allQuality=table();allBands=table();
for i=1:numel(cfg.cases)
    name=cfg.cases{i};
    fprintf('\n[E06-B] %s\n',name);
    path=fullfile(cfg.previousResults,name,'E06_record.mat');
    if ~isfile(path)
        error('E06B:NoRecord', '缺少%s，请先运行E06-A并修改E06B_config.previousResults。',path);
    end
    s=load(path,'Y','Yf','modes','details','t','cfg','info');
    assert(isfield(s,'Y') && isfield(s,'Yf') && isfield(s,'modes'), ...
        'E06-A MAT缺少Y、Yf或modes变量');
    assert(s.info.converged,'已有MVMD分解未收敛，禁止用于改进对照');
    assert(s.details.class_0drop_1mix_2keep(1)==2, ...
        'IMF1并未完整保留，不适合此次冻结的局部改进方案');
    X=s.Y;baseline=s.Yf;imf1=s.modes(:,:,1);t=s.t;
    fs=s.cfg.fs;N=size(X,1);
    [imf1new,shrinkStats]=E06B_joint_haar_denoise(imf1,cfg.detailLevels,cfg.thresholdFactor);
    candidate=baseline+(imf1new-imf1);
    assert(isequal(size(X),size(candidate)) && all(isfinite(candidate(:))), ...
        '输出形状或数值异常');
    parts=regexp(name,'^Z([0-9]+)-([0-9]+)$','tokens','once');
    rollHz=str2double(parts{1});
    switch str2double(parts{2})
        case 1, fm=1;
        case 2, fm=2;
        otherwise,error('E06B:Unsupported','当前仅测试Z3-1和Z3-2');
    end
    freqs=[rollHz,abs(rollHz-fm),rollHz+fm];
    bandNames={'carrier','lower_sideband','upper_sideband'};
    channels={'a','s'};
    for c=1:2
        o=localHarmonicFit(X(:,c),t,freqs);
        b=localHarmonicFit(baseline(:,c),t,freqs);
        q=localHarmonicFit(candidate(:,c),t,freqs);
        oldResidual=sqrt(mean(b.residual.^2));
        newResidual=sqrt(mean(q.residual.^2));
        gain=20*log10(max(oldResidual,eps)/max(newResidual,eps));
        allQuality=[allQuality;table(string(name),string(channels{c}),N,...
            shrinkStats.sigma,shrinkStats.changedPairs,shrinkStats.deltaRMS,...
            oldResidual,newResidual,gain,...
            'VariableNames',{'caseName','channel','samples','estimatedSigmaIMF1',...
            'modifiedHaarPairs','modificationRMS','baselineNonTargetResidualRMS_proxy',...
            'candidateNonTargetResidualRMS_proxy','extraSuppressionDb_proxy'})]; %#ok<AGROW>
        for j=1:3
            r=abs(q.coeff(j))/max(abs(o.coeff(j)),eps)*100;
            bRet=abs(b.coeff(j))/max(abs(o.coeff(j)),eps)*100;
            ph=localWrap(atan2d(imag(q.coeff(j)),real(q.coeff(j)))-...
                atan2d(imag(o.coeff(j)),real(o.coeff(j))));
            allBands=[allBands;table(string(name),string(channels{c}),...
                string(bandNames{j}),freqs(j),bRet,r,ph,...
                'VariableNames',{'caseName','channel','component','frequencyHz',...
                'baselineRetentionPct','candidateRetentionPct','candidatePhaseShiftFromRawDeg'})]; %#ok<AGROW>
        end
    end
    save(fullfile(cfg.outputRoot,[name '_E06B_record.mat']),...
        'X','baseline','candidate','imf1','imf1new','shrinkStats','t','freqs','-v7');
    if i==1
        fig=figure('Visible','off','Color','w','Position',[100 100 1000 690]);
        for c=1:2
            subplot(2,1,c);
            [f,p1]=localSpectrum(X(:,c),fs);
            [~,p2]=localSpectrum(baseline(:,c),fs);
            [~,p3]=localSpectrum(candidate(:,c),fs);
            plot(f,p1,'Color',[.7 .7 .7],'DisplayName','raw');hold on;
            plot(f,p2,'b','DisplayName','frozen MVMD');
            plot(f,p3,'r--','DisplayName','MVMD+joint shrink');
            xlim([0 25]);grid on;
            title(sprintf('%s channel %s',name,channels{c}));ylabel('Amplitude');
            if c==1,legend('Location','northeast');end
            for j=1:3,xline(freqs(j),':','HandleVisibility','off');end
        end
        xlabel('Frequency (Hz)');
        saveas(fig,fullfile(cfg.outputRoot,'E06B_spectrum_preview_Z3-1.png'));close(fig);
    end
end
writetable(allQuality,fullfile(cfg.outputRoot,'E06B_quality.csv'));
writetable(allBands,fullfile(cfg.outputRoot,'E06B_sidebands.csv'));
mask=allBands.component~="carrier";
bandPass=all(allBands.candidateRetentionPct(mask)>=cfg.minRetentionPct & ...
    abs(allBands.candidatePhaseShiftFromRawDeg(mask))<=cfg.maxSidebandPhaseShiftDeg);
noisePass=all(allQuality.extraSuppressionDb_proxy>=cfg.minExtraSuppressionDb);
fid=fopen(fullfile(cfg.outputRoot,'E06B_decision.txt'),'w');
fprintf(fid,'Decision (pre-frozen): sidebands_pass=%d; extra_noise_proxy_pass=%d\n',bandPass,noisePass);
fprintf(fid,'Strong candidate for deeper validation only if both are true.\n');
fprintf(fid,'These measurements are NOT true clean-signal SNR nor attitude RMSE.\n');
fprintf(fid,'DO NOT RETUNE ON Z3-2. If fail, stop this branch.\n');
fclose(fid);
fprintf('\nE06-B结束。动态边带阈值=%d; 额外非目标残差代理阈值=%d\n',bandPass,noisePass);
fprintf('结果在 %s\n',cfg.outputRoot);
end

function h=localHarmonicFit(y,t,freqs)
y=y(:);t=t(:);X=ones(numel(y),1);
for f=freqs
    X=[X,cos(2*pi*f*t),sin(2*pi*f*t)]; %#ok<AGROW>
end
beta=X\y;
h.coeff=complex(beta(2:2:end),beta(3:2:end));
h.residual=y-X*beta;
end
function x=localWrap(x),x=mod(x+180,360)-180;end
function [f,A]=localSpectrum(x,fs)
x=x(:)-mean(x);N=numel(x);half=floor(N/2);
Y=fft(x);A=2*abs(Y(1:half+1))/N;A(1)=A(1)/2;
if rem(N,2)==0,A(end)=A(end)/2;end
f=(0:half)'*fs/N;
end
