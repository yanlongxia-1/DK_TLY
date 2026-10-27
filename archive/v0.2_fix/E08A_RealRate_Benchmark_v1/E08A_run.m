function E08A_run()
%E08A_RUN 真实W类等效有符号角速度的统一滤波对照（非机体全姿态）。
clc;cfg=E08A_config();root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'lib'));
E08A_selftest();
if ~isfolder(cfg.dataRoot)
    error('E08A:DataRoot','缺少原始TXT目录：%s。请先修改E08A_config.m的cfg.dataRoot。',cfg.dataRoot);
end
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
cal=readtable(cfg.calibrationCsv,'TextType','string');
locked=readtable(cfg.lockedParamsCsv,'TextType','string');
requireCols={'operatingMode','rollHz','biasA','biasS','scaleACos','scaleASin',...
    'scaleSCos','scaleSSin','biasSourceFile','recommended'};
assert(all(ismember(requireCols,cal.Properties.VariableNames)),'标定表缺少必需列');
assert(all(ismember({'objective','K','alpha'},locked.Properties.VariableNames)),'锁定参数表缺少字段');
methods=methodList(locked,cfg);
summary=table();statuses=table();
for i=1:numel(cfg.cases)
    filename=char(cfg.cases{i});
    fprintf('\n=== E08-A REAL W: %s ===\n',filename);
    try
        row=runCase(filename,cal,methods,cfg);
        if isempty(summary),summary=row;else,summary=[summary;row];end %#ok<AGROW>
        stat=table(string(filename),"ok", "completed",...
           'VariableNames',{'fileName','status','reason'});
        if isempty(statuses),statuses=stat;else,statuses=[statuses;stat];end %#ok<AGROW>
    catch ME
        msg=sprintf('%s: %s',ME.identifier,ME.message);
        stat=table(string(filename),"failed",string(msg),...
           'VariableNames',{'fileName','status','reason'});
        if isempty(statuses),statuses=stat;else,statuses=[statuses;stat];end %#ok<AGROW>
        fprintf(2,'FAILED %s: %s\n',filename,msg);
        if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
        id=regexprep(filename,'[^A-Za-z0-9]','_');
        fid=fopen(fullfile(cfg.outputDir,['ERROR_' id '.txt']),'w','n','UTF-8');
        if fid>=0,fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));fclose(fid);end
    end
end
writetable(statuses,fullfile(cfg.outputDir,'E08A_status.csv'));
if ~isempty(summary)
    writetable(summary,fullfile(cfg.outputDir,'E08A_summary.csv'));
end
if ~any(strcmp(statuses.status,"ok"))
    fprintf(2,'没有完成真实W文件。先核对dataRoot以及W6-20/W9-20和E6/E9是否齐全。\n');
end
fprintf('E08-A已结束。输出目录：%s\n',cfg.outputDir);
end

function methods=methodList(locked,cfg)
    old=locked(strcmpi(locked.objective,'entropy'),:);
    joint=locked(strcmpi(locked.objective,'joint'),:);
    assert(height(old)==1&&height(joint)==1,'锁定参数必须各有一组entropy与joint');
    methods=struct('name',{},'K',{},'alpha',{},'policy',{});
    methods(1)=struct('name','raw','K',0,'alpha',NaN,'policy','none');
    methods(2)=struct('name','frozen_hard','K',joint.K,'alpha',joint.alpha,'policy','hard');
    methods(3)=struct('name','entropy_hard','K',old.K,'alpha',old.alpha,'policy','hard');
    methods(4)=struct('name','joint_soft','K',joint.K,'alpha',joint.alpha,'policy','soft');
    if cfg.methods.entropyPoolDiagnostic
        methods(5)=struct('name','pool_entropy_hard','K',9,'alpha',450,'policy','hard');
    end
end

function rows=runCase(fileName,cal,methods,cfg)
    tok=regexp(upper(fileName),'^W(\d+)-(\d+)\.TXT$','tokens','once');
    assert(~isempty(tok),'本实验仅支持W滚转+外框恒速工况');
    froll=str2double(tok{1});nominalOuter=str2double(tok{2});
    c=cal(strcmpi(cal.operatingMode,'rolling') & cal.rollHz==froll,:);
    assert(height(c)==1 && logical(c.recommended),'该频率没有唯一且通过验证的固定标定参数');
    coeff=[c.biasA,c.biasS,c.scaleACos,c.scaleASin,c.scaleSCos,c.scaleSSin];
    assert(all(isfinite(coeff)),'标定表含无效值');
    srcW=findUniqueTxt(cfg.dataRoot,fileName);
    srcE=findUniqueTxt(cfg.dataRoot,char(c.biasSourceFile));
    W=read_gyro_txt(srcW);E=read_gyro_txt(srcE);
    assert(W.invalidRows==0 && E.invalidRows==0,'含无效记录行，不能静默改变时间轴');
    assert(W.headerBlockRuleMatches && E.headerBlockRuleMatches,'TXT文件头块计数检查失败');
    N=size(W.data,1);va=evalIndices(N,cfg);
    X=double(W.data(va,[4 7]));tAll=(0:N-1)'/cfg.fs;t=tAll(va);
    outer=rad2deg(unwrap(deg2rad(W.data(:,1))));
    refFit=robust_irls([ones(numel(va),1),t],outer(va),cfg);
    rateTruth=refFit.beta(2); % ONLY evaluation label, NEVER correction input
    inner=rad2deg(unwrap(deg2rad(W.data(:,3))));
    [thetaFixed,phaseConvention]=fixedPhase(inner,tAll,va,froll,cfg);
    thetaEncoder=inner(va); % offline diagnostic ONLY, not independent roll estimation
    [eStart,eLast]=trimE(size(E.data,1),cfg);
    fa=robust_harmonic_fit(E.data(eStart:eLast,3),E.data(eStart:eLast,4),cfg);
    fs=robust_harmonic_fit(E.data(eStart:eLast,3),E.data(eStart:eLast,7),cfg);
    HE=[fa.cosCoefficient fa.sinCoefficient;fs.cosCoefficient fs.sinCoefficient];
    bE=[c.biasA c.biasS];
    scale=[c.scaleACos c.scaleASin;c.scaleSCos c.scaleSSin];
    [~,stem]=fileparts(fileName);out=fullfile(cfg.outputDir,stem);
    if ~isfolder(out),mkdir(out);end
    % 每个唯一K/alpha联合MVMD仅运行一次，两个重构规则复用同一份模态。
    data=struct();data.raw=X;solver=table();cache=struct('key',{},'hard',{},'soft',{});
    for mi=2:numel(methods)
        mm=methods(mi);key=sprintf('K%d_A%.9f',mm.K,mm.alpha);
        idx=find(strcmp({cache.key},key),1);
        if isempty(idx)
            fprintf('  MVMD联合分解: K=%d alpha=%.6f ...\n',mm.K,mm.alpha);
            tic;
            [modes,centerHz,info]=gyro_mvmd(X,cfg.fs,mm.K,mm.alpha,...
                cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
            elapsed=toc;
            valid=logical(info.converged && ...
                info.AC_reconRelative<=cfg.mvmd.maxAcReconstructionLoss);
            solverRow=table(string(stem),mm.K,mm.alpha,...
                logical(info.converged),info.iterations,info.AC_reconRelative,elapsed,valid,...
                'VariableNames',{'fileName','K','alpha','converged','iterations',...
                'ACreconRelative','elapsedSeconds','valid'});
            if isempty(solver),solver=solverRow;else,solver=[solver;solverRow];end %#ok<AGROW>
            if valid
                localCfg=cfg;localCfg.rollHz=froll;
                [hard,detailsH]=gyro_reconstruct(modes,centerHz,X,localCfg,'paper');
                [soft,detailsS]=E07D_soft_reconstruct(modes,centerHz,X,localCfg);
                detailsH.source= repmat("hard",height(detailsH),1);
                detailsS.source= repmat("soft",height(detailsS),1);
                detailsH.K=repmat(mm.K,height(detailsH),1);
                detailsH.alpha=repmat(mm.alpha,height(detailsH),1);
                detailsS.K=repmat(mm.K,height(detailsS),1);
                detailsS.alpha=repmat(mm.alpha,height(detailsS),1);
                % 存单独表避免hard/soft不同列结构产生拼接异常
                writetable(detailsH,fullfile(out,sprintf('modes_K%d_A%g_hard.csv',mm.K,mm.alpha)));
                writetable(detailsS,fullfile(out,sprintf('modes_K%d_A%g_soft.csv',mm.K,mm.alpha)));
            else
                hard=[];soft=[];
                fprintf(2,'  该参数未收敛或重构不合格: K=%d alpha=%g\n',mm.K,mm.alpha);
            end
            cache(end+1)=struct('key',key,'hard',hard,'soft',soft); %#ok<AGROW>
            idx=numel(cache);
        end
        if strcmp(mm.policy,'hard'),yy=cache(idx).hard;else,yy=cache(idx).soft;end
        if ~isempty(yy),data.(mm.name)=yy;end
    end
    writetable(solver,fullfile(out,'E08A_solver.csv'));
    methodsAvailable=fieldnames(data);metrics=table();T=table(t,thetaFixed,thetaEncoder,...
        repmat(rateTruth,numel(t),1),'VariableNames',{'time_s','phase_fixed_deg',...
        'phase_encoder_diagnostic_deg','outer_encoder_reference_deg_s'});
    for mi=1:numel(methodsAvailable)
        method=methodsAvailable{mi};yy=data.(method);
        [rateFixed,den]=E08A_estimate_outer(yy,thetaFixed,HE,bE,scale);
        [rateEnc,~]=E08A_estimate_outer(yy,thetaEncoder,HE,bE,scale);
        fit=carrierMetrics(X,yy,t,froll);
        err=rateFixed-rateTruth;errEnc=rateEnc-rateTruth;
        rowsNow=table(string(stem),string(method),froll,nominalOuter,rateTruth,...
            mean(rateFixed),mean(err),std(err),mean(abs(err)),sqrt(mean(err.^2)),...
            sqrt(mean(errEnc.^2)),min(den),fit.keepA,fit.keepS,...
            fit.phaseShiftA,fit.phaseShiftS,fit.relativePhaseShift,...
            fit.snrProxyGainA,fit.snrProxyGainS,numel(t),...
            'VariableNames',{'fileName','method','rollHz','nominalOuterDegS',...
            'referenceOuterDegS','meanEstimatedDegS','meanErrorDegS','errorStdDegS',...
            'MAE_deg_s','RMSE_deg_s','RMSE_encoderPhaseDiagnostic_deg_s',...
            'minResponsePower','carrierRetentionA_pct','carrierRetentionS_pct',...
            'carrierPhaseShiftA_deg','carrierPhaseShiftS_deg',...
            'relativePhaseShift_deg','snrProxyGainA_db','snrProxyGainS_db',...
            'evaluationSamples'});
        if isempty(metrics),metrics=rowsNow;else,metrics=[metrics;rowsNow];end %#ok<AGROW>
        T.([method '_fixed_deg_s'])=rateFixed;
        T.([method '_encoder_diagnostic_deg_s'])=rateEnc;
    end
    rows=metrics;
    writetable(metrics,fullfile(out,'E08A_rate_metrics.csv'));
    writetable(T,fullfile(out,'E08A_rate_series.csv'));
    save(fullfile(out,'E08A_record.mat'),'data','metrics','T','solver','HE','bE','scale',...
        'cfg','c','phaseConvention','-v7');
    plotRate(T,methodsAvailable,stem,out,cfg);
    fid=fopen(fullfile(out,'E08A_interpretation.txt'),'w','n','UTF-8');
    if fid>=0
        fprintf(fid,'真实W类恒速外框等效角速度（deg/s），不是一般机体Omega_x/Omega_y或姿态角RMSE。\n');
        fprintf(fid,'主解调用已知滚转频率+初始角和符号；相位初始化方式: %s。\n',phaseConvention);
        fprintf(fid,'所有方案共用E类频率组固定偏置、E类零输入谐波、复响应标定。没有用评价参考来调参。\n');
        fprintf(fid,'此处没有补偿W类工况有效截距差值，它可能造成滚转频率上的假角速度。\n');
        fprintf(fid,'完整分块MVMD为离线非因果实现；W6-20可能已参与之前算法开发，不当作全新盲测。\n');
        fprintf(fid,'pool_entropy_hard是E07D共同候选池事后熵最小诊断，不是独立FVIM搜索胜者。\n');
        fclose(fid);
    end
    fprintf('  %s参考 %+0.4f deg/s；完成%u种算法输出。\n',stem,rateTruth,height(metrics));
end

function [thetaFixed,convention]=fixedPhase(inner,tAll,va,froll,cfg)
    if strcmpi(cfg.phase.mode,'encoder_initial')
        phi0=inner(1);
        signRoll=sign(median(diff(inner)));
        assert(signRoll~=0,'内框角度无法确定滚转方向');
        convention='bench_encoder_initial_phase_and_direction';
    elseif strcmpi(cfg.phase.mode,'manual')
        phi0=cfg.phase.manualPhi0Deg;signRoll=cfg.phase.manualRollSign;
        assert(isfinite(phi0)&&ismember(signRoll,[-1,1]),'manual模式下设置phi0和+/-方向');
        convention='user_precalibrated_initial_phase_and_direction';
    else
        error('不支持的初始相位来源');
    end
    thetaFixed=phi0+signRoll*360*froll*tAll(va);
end
function [a,b]=trimE(N,cfg)
    a=1+round(cfg.EtrimStartSec*cfg.fs);b=N-round(cfg.EtrimEndSec*cfg.fs);
    if b-a+1<round(cfg.minimumEvalSec*cfg.fs),a=1;b=N;end
end
function idx=evalIndices(N,cfg)
    idx=(max(1,round(cfg.evalFraction(1)*N)):...
        min(N,round(cfg.evalFraction(2)*N)))';
    assert(numel(idx)>=round(cfg.minimumEvalSec*cfg.fs),...
        '评价数据不足最短时间，请选择记录更长的W文件');
end
function fp=findUniqueTxt(root,file)
    d=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'))];
    ids=arrayfun(@(x)strcmpi(x.name,file),d);d=d(ids);
    if ~isempty(d)
        p=fullfile({d.folder},{d.name});pl=cellfun(@lower,p,'UniformOutput',false);
        [~,ix]=unique(pl,'stable');d=d(ix);
    end
    assert(numel(d)==1,'需要在dataRoot下找到唯一的%s（目前找到%d个）',file,numel(d));
    fp=fullfile(d.folder,d.name);
end
function f=carrierMetrics(X,Y,t,rollHz)
    f0=fitTwo(X,t,rollHz);f1=fitTwo(Y,t,rollHz);
    f.keepA=100*f1.amp(1)/max(f0.amp(1),eps);
    f.keepS=100*f1.amp(2)/max(f0.amp(2),eps);
    da=wrap180(f1.phase(1)-f0.phase(1));
    ds=wrap180(f1.phase(2)-f0.phase(2));
    f.phaseShiftA=da;f.phaseShiftS=ds;
    f.relativePhaseShift=wrap180(ds-da);
    f.snrProxyGainA=f1.proxy(1)-f0.proxy(1);
    f.snrProxyGainS=f1.proxy(2)-f0.proxy(2);
end
function f=fitTwo(X,t,rollHz)
    phi=2*pi*rollHz*t;
    B=[ones(size(t)),cos(phi),sin(phi)];
    q=B\X;
    model=B*q;resid=X-model;
    f.amp=hypot(q(2,:),q(3,:));
    f.phase=atan2d(q(3,:),q(2,:));
    ac=model-repmat(q(1,:),size(model,1),1);
    f.proxy=10*log10(sum(ac.^2,1)./max(sum(resid.^2,1),eps));
end
function x=wrap180(x)
    x=mod(x+180,360)-180;
end
function plotRate(T,methods,stem,out,cfg)
    fig=figure('Visible',cfg.figureVisible,'Color','w','Position',[100 100 1250 660]);
    plot(T.time_s,T.outer_encoder_reference_deg_s,'k--','LineWidth',1.4,'DisplayName','outer encoder');hold on;
    colors=lines(numel(methods));
    for i=1:numel(methods)
        nm=methods{i};plot(T.time_s,T.([nm '_fixed_deg_s']),'Color',colors(i,:),...
            'DisplayName',strrep(nm,'_',' '));
    end
    xlabel('time (s)');ylabel('signed outer equivalent rate (deg/s)');
    title([stem ' | constant-roll-frequency demodulation']);grid on;
    legend('Location','best');
    saveas(fig,fullfile(out,'E08A_rate_compare.png'));close(fig);
end
