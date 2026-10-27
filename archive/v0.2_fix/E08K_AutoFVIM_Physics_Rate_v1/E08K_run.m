function E08K_run()
%E08K_RUN 自动FVIM风格寻优 vs 固定强基准，训练不接触转台真值。
% 不使用每频率K/alpha锁定表。两优化目标各4个提案、共享4个初始候选。
% 原论文FVIM算法细节未完全核实，当前搜索器是四领先者工程复现，不冒称作者源码。
clc;root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
cfg=E08K_config();E08K_selftest();
assert(isfolder(cfg.dataRoot),'E08K:DataMissing','修改E08K_config.m中cfg.dataRoot，当前不存在: %s',cfg.dataRoot);
assert(size(cfg.search.initial,1)==4 && cfg.search.proposalsPerObjective>=1);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
save(fullfile(cfg.outputDir,'E08K_frozen_config.mat'),'cfg');
rng(cfg.search.seed,'twister');
allCandidates=table();allSolvers=table();allTraces=table();allLocked=table();allEvaluation=table();allEvalSolvers=table();allStatus=table();
for fi=1:numel(cfg.frequenciesHz)
    f=cfg.frequenciesHz(fi);
    fprintf('\n=== E08-K %g Hz | 自动寻优: 不查询K-alpha频率表 ===\n',f);
    try
        dTrain={loadData(sprintf('W%d-10.TXT',f),f,'train',cfg),...
                loadData(sprintf('Z%d-1.TXT',f),f,'train',cfg)};
    catch ME
        warning('E08K:TrainingFile','%gHz训练文件不可用: %s',f,ME.message);
        allStatus=addRow(allStatus,table(f,string('training'),string('failed'),string(ME.message),...
            'VariableNames',{'rollHz','stage','status','message'}));persist();continue;
    end
    cache=struct('key',{},'K',{},'alpha',{},'oldFit',{},'newFit',{},'oldFeasible',{},'newFeasible',{});
    oldPool=struct('K',{},'alpha',{},'score',{});newPool=oldPool;
    init=cfg.search.initial;
    for ii=1:4
        c=evaluate(init(ii,1),init(ii,2),'shared_seed');
        oldPool(end+1)=struct('K',c.K,'alpha',c.alpha,'score',c.oldFit); %#ok<AGROW>
        newPool(end+1)=struct('K',c.K,'alpha',c.alpha,'score',c.newFit); %#ok<AGROW>
    end
    for g=1:cfg.search.proposalsPerObjective
        [k,a]=propose(oldPool,g,cfg,cache);
        old=evaluate(k,a,'entropy_proposal');
        oldPool(end+1)=struct('K',old.K,'alpha',old.alpha,'score',old.oldFit); %#ok<AGROW>
        [k,a]=propose(newPool,g,cfg,cache);
        joint=evaluate(k,a,'joint_proposal');
        newPool(end+1)=struct('K',joint.K,'alpha',joint.alpha,'score',joint.newFit); %#ok<AGROW>
        [bestO,validO]=bestOf(oldPool);[bestN,validN]=bestOf(newPool);
        trace=table(f,g,numel(cache),validO,validN,bestO.K,bestO.alpha,bestO.score,...
            bestN.K,bestN.alpha,bestN.score,...
            'VariableNames',{'rollHz','generation','uniqueCandidates','oldFeasible','newFeasible',...
            'oldBestK','oldBestAlpha','oldBestFitness','jointBestK','jointBestAlpha','jointBestFitness'});
        allTraces=addRow(allTraces,trace);persist();
    end
    % IMPORTANT: select each objective from SAME UNION of ALL evaluated candidates.
    % The proposer populations are independent, but assessment is paired on shared cache.
    [old,oldOK]=bestFromCommonPool(cache,'oldFit');
    [joint,newOK]=bestFromCommonPool(cache,'newFit');
    methods=struct('name',{},'K',{},'alpha',{},'enabled',{});
    methods(1)=struct('name','raw','K',0,'alpha',NaN,'enabled',true);
    methods(2)=struct('name','fixed_baseline','K',8,'alpha',154.157959862,'enabled',true);
    methods(3)=struct('name','developer_K5','K',5,'alpha',154.157959862,'enabled',true);
    methods(4)=struct('name','entropy_auto','K',old.K,'alpha',old.alpha,'enabled',oldOK);
    methods(5)=struct('name','physics_auto','K',joint.K,'alpha',joint.alpha,'enabled',newOK);
    for mm=2:5
        m=methods(mm);
        allLocked=addRow(allLocked,table(f,string(m.name),m.K,m.alpha,m.enabled,...
            'VariableNames',{'rollHz','method','K','alpha','trainingFeasible'}));
    end
    persist();
    zmod=min(cfg.testZmodulation,f-1);
    names={sprintf('W%d-%d.TXT',f,cfg.testWspeed),sprintf('Z%d-%d.TXT',f,zmod)};
    if cfg.runDevelopmentDiagnostic
        names{end+1}=sprintf('W%d-20.TXT',f); %#ok<AGROW>
    end
    for testid=1:numel(names)
        fname=names{testid};isDevelopment=testid>2;
        try
            d=loadData(fname,f,'test',cfg);
            if isDevelopment,d.validationRole="development_diagnostic";else,d.validationRole="locked_external_check";end
            [rows,sols]=evaluateHoldout(d,methods,cfg);
            allEvaluation=addRows(allEvaluation,rows);
            allEvalSolvers=addRows(allEvalSolvers,sols);
            allStatus=addRow(allStatus,table(f,string(fname),string('evaluated'),string(''),...
                'VariableNames',{'rollHz','stage','status','message'}));
            fprintf('  已锁定参数后检查: %s，%d条方法记录。\n',fname,height(rows));
        catch ME
            warning('E08K:Holdout','%s 未完成检查：%s',fname,ME.message);
            allStatus=addRow(allStatus,table(f,string(fname),string('failed'),string(ME.message),...
                'VariableNames',{'rollHz','stage','status','message'}));
        end
        persist();
    end
end
writeDecision(cfg,allLocked,allEvaluation,allStatus);
fprintf('\nE08-K结束。结果：%s\n',cfg.outputDir);

    function c=evaluate(K,alpha,origin)
        K=round(K);alpha=round(alpha*1e4)/1e4;
        key=sprintf('K%d_A%.4f',K,alpha);
        jj=find(strcmp({cache.key},key),1);
        if ~isempty(jj),c=cache(jj);return;end
        assert(numel(cache)<cfg.search.maxUniqueCandidates,'E08K:Budget','超出每频率MVMD候选数量上限');
        o=nan(2,1);n=nan(2,1);solverValid=false(2,1);phyValid=false(2,1);
        for ti=1:2
            r=solveOne(dTrain{ti},K,alpha,cfg);
            solverValid(ti)=r.solverValid;
            phyValid(ti)=r.physicalPass;
            o(ti)=r.entropyFit;n(ti)=r.jointFit;
            allSolvers=addRow(allSolvers,solverRow(r,string(origin),'TRAIN'));
            persist();
        end
        oldOK=all(solverValid) && all(isfinite(o));
        newOK=oldOK && all(phyValid) && all(isfinite(n));
        oldFit=Inf;newFit=Inf;
        if oldOK,oldFit=mean(o);end
        if newOK,newFit=mean(n);end
        c=struct('key',key,'K',K,'alpha',alpha,'oldFit',oldFit,'newFit',newFit,...
            'oldFeasible',oldOK,'newFeasible',newOK);
        cache(end+1)=c; %#ok<AGROW>
        allCandidates=addRow(allCandidates,table(f,K,alpha,string(origin),oldOK,newOK,oldFit,newFit,...
            'VariableNames',{'rollHz','K','alpha','source','entropyFeasible','jointFeasible','entropyFitness','jointFitness'}));
        persist();
        fprintf('[%g Hz #%d/%d] K=%d alpha=%.4f | entropy=%g joint=%g\n',f,numel(cache),...
            cfg.search.maxUniqueCandidates,K,alpha,oldFit,newFit);
    end
    function persist()
        writeTable(allCandidates,'E08K_train_candidates.csv');
        writeTable(allSolvers,'E08K_train_solver.csv');
        writeTable(allTraces,'E08K_search_trace.csv');
        writeTable(allLocked,'E08K_locked_params.csv');
        writeTable(allEvaluation,'E08K_test_metrics.csv');
        writeTable(allEvalSolvers,'E08K_test_solver.csv');
        writeTable(allStatus,'E08K_status.csv');
    end
    function writeTable(T,name)
        if ~isempty(T),writetable(T,fullfile(cfg.outputDir,name));end
    end
end

function d=loadData(filename,f,role,cfg)
path=findUnique(cfg.dataRoot,filename);r=read_gyro_txt(path);
assert(r.invalidRows==0 && r.headerBlockRuleMatches,'TXT异常记录或头信息错误：%s',filename);
N=size(r.data,1);assert(N>cfg.fs,'原始文件过短：%s',filename);
% Frozen E08-J last-window rule applies BOTH train and test; no per-file cherry-picking.
last=min(N,floor(cfg.endFraction*N));
first=max(1+round(cfg.skipStartSec*cfg.fs),last-cfg.lastSeconds*cfg.fs+1);
ix=(first:last)';
assert(numel(ix)>=cfg.minTestSeconds*cfg.fs,'有效窗口不足%.1f秒：%s',cfg.minTestSeconds,filename);
d.X=double(r.data(ix,[4 7]));d.t=(0:numel(ix)-1)'/cfg.fs;
d.fileName=filename;d.froll=f;d.index=ix;d.all=double(r.data);d.absoluteTime=(ix-1)/cfg.fs;
d.kind=upper(filename(1));
if d.kind=='Z'
    tok=regexp(upper(filename),'^Z\d+-(\d+)\.TXT$','tokens','once');
    assert(~isempty(tok),'Z文件命名不符合规则');fm=str2double(tok{1});
    d.freqs=[f-fm f f+fm];d.side=[f-fm f+fm];
else
    d.freqs=f;d.side=f;
end
assert(all(d.freqs>0),'调制频率应小于滚转频率，当前文件：%s',filename);
end

function path=findUnique(root,filename)
% Windows不区分大小写；扫描全部TXT并用精确basename检查重复。
a=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'));
   dir(fullfile(root,'*.TXT'));dir(fullfile(root,'*.txt'))];
ids=arrayfun(@(z) strcmpi(z.name,filename),a);a=a(ids);
if ~isempty(a)
    paths=fullfile({a.folder},{a.name});
    [~,u]=unique(lower(string(paths)),'stable');a=a(u);
end
assert(numel(a)==1,'必须找到唯一文件 %s，目前匹配%d个',filename,numel(a));
path=fullfile(a.folder,a.name);
end

function r=solveOne(d,K,alpha,cfg)
% Strict original MVMD solver (eta=1), FIXED E08-J amplitude/phase criterion.
% No encoder truth or calibrated angular-rate value enters either fitness.
tStart=tic;r.fileName=string(d.fileName);r.rollHz=d.froll;r.K=K;r.alpha=alpha;
r.solverValid=false;r.physicalPass=false;r.entropyFit=Inf;r.jointFit=Inf;
r.iterations=NaN;r.delta=NaN;r.centerDeltaHz=NaN;r.acLoss=NaN;r.converged=false;
r.suppA=NaN;r.suppS=NaN;r.retMin=NaN;r.retMax=NaN;r.phaseMax=NaN;
r.Y=[];r.errorMessage="";
try
    [modes,cent,info]=E08G_mvmd_damped(d.X,cfg.fs,K,alpha,...
        cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau,d.froll,1,[]);
    r.iterations=info.iterations;r.delta=info.delta;r.centerDeltaHz=info.centerDeltaHz;
    r.acLoss=info.ACreconRelative;r.converged=logical(info.converged);
    r.solverValid=logical(info.converged && info.ACreconRelative<=cfg.mvmd.maxAcReconstructionLoss);
    if r.solverValid
        local=cfg;local.rollHz=d.froll;
        [Y,~]=gyro_reconstruct(modes,cent,d.X,local,'paper');
        [q,band]=E08J_physics_metrics(d.X,Y,d.t,d.freqs,local);
        r.entropyFit=E07D_entropy(modes);
        % Both channels must improve in non-target residual proxy; noise proxy is NOT truth SNR.
        r.suppA=band.proxyDb(find(band.channel==1,1));
        r.suppS=band.proxyDb(find(band.channel==2,1));
        r.retMin=q.minimumStrongRetention;r.retMax=q.maximumStrongRetention;
        r.phaseMax=max(q.maximumStrongPhaseDeg,q.maximumStrongRelativePhaseDeg);
        r.physicalPass=logical(q.newPass && r.suppA>=0 && r.suppS>=0);
        if r.physicalPass
            r.jointFit=-(q.minimumProxyDb+0.1*q.meanProxyDb-0.001*K);
        end
        r.Y=Y;
    end
catch ME
    r.errorMessage=string(ME.message);
    warning('E08K:SolverError','%s K%d alpha%.4f: %s',d.fileName,K,alpha,ME.message);
end
r.seconds=toc(tStart);
end

function row=solverRow(r,role,phase)
row=table(r.fileName,r.rollHz,string(role),string(phase),r.K,r.alpha,...
    r.converged,r.iterations,r.delta,r.centerDeltaHz,r.acLoss,...
    r.solverValid,r.physicalPass,r.seconds,...
    'VariableNames',{'fileName','rollHz','role','stage','K','alpha','converged',...
    'iterations','delta','centerDeltaHz','ACreconRelative','solverValid','physicalPass','seconds'});
end

function [K,alpha]=propose(pop,gen,cfg,cache)
% FVIM风格四领先者，两个目标各自维护自己的搜索历史/种群。
lo=[cfg.search.Kmin log10(cfg.search.alphaMin)];
hi=[cfg.search.Kmax log10(cfg.search.alphaMax)];
scores=[pop.score];[~,ix]=sort(scores,'ascend');
ll=pop(ix(1:min(4,numel(ix))));
u=[[ll.K]' log10([ll.alpha]')];
a=2*(1-gen/(cfg.search.proposalsPerObjective+1));
for attempt=1:80
    if all(~isfinite(scores)) || attempt>25
        p=lo+rand(1,2).*(hi-lo);
    else
        center=mean(u,1);mv=zeros(size(u));
        for j=1:size(u,1)
            A=a*(2*rand(1,2)-1);C=2*rand(1,2);
            mv(j,:)=u(j,:)-A.*abs(C.*u(j,:)-center);
        end
        p=mean(mv,1)+.06*randn(1,2).*(hi-lo);
    end
    p=max(lo,min(hi,p));K=round(p(1));alpha=round(10^p(2)*1e4)/1e4;
    if ~any([cache.K]==K & abs([cache.alpha]-alpha)<1e-5),return;end
end
error('E08K:CandidatesExhausted','无法产生未被缓存的新候选');
end

function [b,valid]=bestOf(pop)
scores=[pop.score];[v,i]=min(scores);valid=isfinite(v);
b=pop(i);if ~valid,b.K=NaN;b.alpha=NaN;b.score=Inf;end
end

function [b,valid]=bestFromCommonPool(cache,scoreField)
% All objectives are re-scored on the same explored union; no best-of-own-path bias.
assert(~isempty(cache),'Empty candidate pool');
score=[cache.(scoreField)];[v,i]=min(score);valid=isfinite(v);
b=struct('K',cache(i).K,'alpha',cache(i).alpha,'score',v);
if ~valid,b.K=NaN;b.alpha=NaN;b.score=Inf;end
end

function T=addRow(T,row)
if isempty(T),T=row;else,T=[T;row];end
end
function T=addRows(T,rows)
if isempty(rows),return;end
if isempty(T),T=rows;else,T=[T;rows];end
end

function [rows,solvers]=evaluateHoldout(d,methods,cfg)
% 参数已在训练阶段锁定：这里严禁更换K/alpha、评分权重或调整标定参数。
rows=table();solvers=table();
cal=readtable(cfg.calibrationCsv,'TextType','string');
C=cal(strcmpi(cal.operatingMode,'rolling') & cal.rollHz==d.froll,:);
assert(height(C)==1 && logical(C.recommended),'该频率没有唯一的可用预标定记录');
param=[C.biasA C.biasS C.scaleACos C.scaleASin C.scaleSCos C.scaleSSin];
assert(all(isfinite(param)),'频率标定参数缺失');
% 仅W类有符号外框等效角速度可以使用外框编码器求RMSE；Z只评物理代理。
truth=NaN;HE=[];bE=[];scale=[];phase=[];
if d.kind=='W'
    assert(strcmpi(cfg.phase.mode,'encoder_initial'),'当前角速度评估只支持encoder_initial相位初始化');
    outer=rad2deg(unwrap(deg2rad(d.all(:,1))));
    fit=robust_irls([ones(numel(d.index),1) d.absoluteTime],outer(d.index),cfg);
    truth=fit.beta(2);
    inner=rad2deg(unwrap(deg2rad(d.all(:,3))));
    sgn=sign(median(diff(inner)));
    assert(abs(sgn)==1,'内框滚转方向无法确定');
    phase=inner(1)+sgn*360*d.froll*d.absoluteTime;
    E=read_gyro_txt(findUnique(cfg.dataRoot,char(C.biasSourceFile)));
    assert(E.invalidRows==0 && E.headerBlockRuleMatches,'E类标定原始TXT不合格');
    ix=(1+round(cfg.EtrimStartSec*cfg.fs)):(size(E.data,1)-round(cfg.EtrimEndSec*cfg.fs));
    assert(numel(ix)>cfg.fs,'E类标定记录过短');
    fa=robust_harmonic_fit(double(E.data(ix,3)),double(E.data(ix,4)),cfg);
    fs=robust_harmonic_fit(double(E.data(ix,3)),double(E.data(ix,7)),cfg);
    HE=[fa.cosCoefficient fa.sinCoefficient;fs.cosCoefficient fs.sinCoefficient];
    bE=[C.biasA C.biasS];
    scale=[C.scaleACos C.scaleASin;C.scaleSCos C.scaleSSin];
end
cached=struct('key',{},'result',{});
T=table(d.absoluteTime,'VariableNames',{'time_s'});
if isfinite(truth),T.reference_outer_deg_s=repmat(truth,numel(d.t),1);end
for i=1:numel(methods)
    method=methods(i);
    if ~method.enabled
        rows=addRow(rows,evalRow(d,method,false,false,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN));
        continue;
    end
    if method.K==0
        r=struct('Y',d.X,'solverValid',true,'physicalPass',true,'suppA',0,'suppS',0,...
            'retMin',100,'retMax',100,'phaseMax',0,'converged',true,'iterations',0,...
            'delta',0,'centerDeltaHz',0,'acLoss',0,'seconds',0,'fileName',string(d.fileName),...
            'rollHz',d.froll,'K',0,'alpha',NaN);
    else
        key=sprintf('K%d_A%.4f',method.K,round(method.alpha*1e4)/1e4);
        ci=find(strcmp({cached.key},key),1);
        if isempty(ci)
            r=solveOne(d,method.K,method.alpha,cfg);
            cached(end+1)=struct('key',key,'result',r); %#ok<AGROW>
            solvers=addRow(solvers,solverRow(r,string(method.name),'TEST'));
        else
            r=cached(ci).result;
        end
    end
    mae=NaN;rmse=NaN;errMean=NaN;errStd=NaN;meanRate=NaN;
    if r.solverValid && ~isempty(r.Y) && isfinite(truth)
        [rate,~]=E08A_estimate_outer(r.Y,phase,HE,bE,scale);
        err=rate-truth;
        mae=mean(abs(err));rmse=sqrt(mean(err.^2));errMean=mean(err);
        errStd=std(err);meanRate=mean(rate);
        T.(char(method.name))=rate;
    end
    rows=addRow(rows,evalRow(d,method,r.solverValid,r.physicalPass,r.suppA,r.suppS,...
        r.retMin,r.retMax,r.phaseMax,truth,meanRate,mae,rmse,errMean,errStd));
end
% 单个测试文件分别保存曲线，后续可计算paired错误差。Z类只有频域代理，无角速度真值。
stem=regexprep(d.fileName,'[^A-Za-z0-9]','_');
writetable(T,fullfile(cfg.outputDir,['E08K_series_' char(stem) '.csv']));
if cfg.makeFigures && isfinite(truth)
    try
        fig=figure('Visible','off','Color','w');hold on;
        plot(T.time_s,T.reference_outer_deg_s,'k--','LineWidth',1.2);
        legendEntries={'encoder reference'};
        for i=1:numel(methods)
            n=char(methods(i).name);
            if ismember(n,T.Properties.VariableNames)
                plot(T.time_s,T.(n));legendEntries{end+1}=n; %#ok<AGROW>
            end
        end
        grid on;xlabel('Time(s)');ylabel('Equivalent signed outer rate (deg/s)');
        title(sprintf('%s: held-out fixed-calibration estimate',d.fileName));
        legend(legendEntries,'Interpreter','none','Location','best');
        saveas(fig,fullfile(cfg.outputDir,['E08K_rate_' char(stem) '.png']));close(fig);
    catch ME
        warning('E08K:Plot','跳过图像：%s',ME.message);
        if exist('fig','var') && isgraphics(fig),close(fig);end
    end
end
end

function row=evalRow(d,m,solverValid,physicalPass,suppA,suppS,retMin,retMax,phaseMax,truth,meanRate,mae,rmse,errMean,errStd)
% 不把SNR代理写作真实SNR；Z类没有外框恒速RMSE，统一NaN。
row=table(string(d.fileName),string(d.validationRole),d.froll,string(m.name),m.K,m.alpha,...
    logical(solverValid),logical(physicalPass),suppA,suppS,retMin,retMax,phaseMax,...
    truth,meanRate,mae,rmse,errMean,errStd,numel(d.t),...
    'VariableNames',{'fileName','validationRole','rollHz','method','K','alpha','solverValid','fidelityPass',...
    'noiseProxyA_dB','noiseProxyS_dB','minRetentionPct','maxRetentionPct','maxPhaseDeg',...
    'encoderReferenceDegS','meanRateDegS','MAE_deg_s','RMSE_deg_s','meanErrorDegS',...
    'stdErrorDegS','evaluationSamples'});
end

function writeDecision(cfg,locked,checks,status)
fid=fopen(fullfile(cfg.outputDir,'E08K_decision.txt'),'w','n','UTF-8');
assert(fid>=0,'无法写入总结文件');cl=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'E08-K：FVIM风格自动双参数搜索（非原论文作者精确源码）。\n');
fprintf(fid,'算法只读取训练双通道和已知频率；测试文件/编码器角速度不参与FVIM评分。\n');
fprintf(fid,'初始4候选共享，旧熵目标/新联合目标各额外%d候选；最多%d组候选/频率。\n',...
 cfg.search.proposalsPerObjective,cfg.search.maxUniqueCandidates);
fprintf(fid,'对照：固定K=8 alpha=154.157959862、开发K5/154.157959862、旧熵目标自动选参、新物理目标自动选参、原始信号。\n');
fprintf(fid,'相同硬重构，比较时必须分别统计成功率及只对同文件有效方法的成对RMSE。\n');
fprintf(fid,'W类外框RMSE使用同频率固定校准表与E类谐波基线；Z类当前仅有滤波代理指标。\n');
fprintf(fid,'噪声是非目标残差单频/边带拟合代理，并非有干净真值的SNR。\n');
fprintf(fid,'W9-20等曾用于开发，不可称盲测；W9-60/Z9-3的历史使用情况亦须审计。\n');
if ~isempty(locked)
    for i=1:height(locked)
        fprintf(fid,'%gHz %s K=%g alpha=%.5g trainingFeasible=%d\n',...
        locked.rollHz(i),char(locked.method(i)),locked.K(i),locked.alpha(i),locked.trainingFeasible(i));
    end
end
if ~isempty(checks)
    for i=1:height(checks)
        fprintf(fid,'TEST %s %s solverValid=%d fidelityPass=%d RMSEdegS=%.7g proxyA=%.4f proxyS=%.4f\n',...
            char(checks.fileName(i)),char(checks.method(i)),checks.solverValid(i),checks.fidelityPass(i),...
            checks.RMSE_deg_s(i),checks.noiseProxyA_dB(i),checks.noiseProxyS_dB(i));
    end
end
if ~isempty(status)
    for i=1:height(status)
        if status.status(i)=="failed"
            fprintf(fid,'MISSING/FAILED %s: %s\n',char(status.stage(i)),char(status.message(i)));
        end
    end
end
fprintf(fid,'不能从有限训练/测试证明全局最优或在线实时性；此处是分段离线自动寻优。\n');
end
