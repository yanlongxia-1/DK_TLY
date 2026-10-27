function E08C_run()
%E08C_RUN 3/9Hz W/Z真实工况的条件化MVMD + FVIM风格联合目标搜索可行性。
% 提醒：FVIM风格四领先者优化是工程复现，非论文作者完整FVIM精确源码。
clc;cfg=E08C_config();root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'lib'));
E08C_selftest();
assert(isfolder(cfg.dataRoot),'原始数据目录不存在：%s',cfg.dataRoot);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
rng(cfg.search.seed,'twister');
solver=table();train=table();locked=table();checks=table();diagTbl=table();
for ii=1:numel(cfg.frequenciesHz)
    f=cfg.frequenciesHz(ii);
    fprintf('\n==== E08-C | 已知滚转 %.0f Hz；训练W/Z，不看测试结果调参 ====\n',f);
    trainFiles={sprintf('W%d%s',f,cfg.trainWsuffix),sprintf('Z%d%s',f,cfg.trainZsuffix)};
    checkFiles={sprintf('W%d%s',f,cfg.checkWsuffix),sprintf('Z%d%s',f,cfg.checkZsuffix)};
    tr=cell(2,1);
    for j=1:2, tr{j}=loadCase(trainFiles{j},f,j==2,cfg);end
    % 第0关：仅比较同一参数的不同初始中心，不改变收敛判据。
    for initFlag=0:1
        r=solve(tr{1},8,154.157959862,logical(initFlag),'init_diagnostic',cfg);
        solver=appendRow(solver,r.solver);persist();
        dd=table(f,string(trainFiles{1}),logical(initFlag),r.solver.valid,r.solver.iterations,...
            r.solver.delta,r.solver.centerDeltaHz,r.solver.ACreconRelative,...
            'VariableNames',{'rollHz','fileName','frequencyInitialized','valid',...
            'iterations','delta','centerDeltaHz','ACreconRelative'});
        diagTbl=appendRow(diagTbl,dd);persist();
    end
    % 仅训练集调 K/alpha；候选的两个文件都要收敛并且满足幅相约束。
    cache=struct('K',{},'alpha',{},'fitness',{},'noiseA',{},'noiseS',{},'minRetention',{},'maxPhase',{},'valid',{});
    initial=cfg.search.initial;
    for q=1:size(initial,1)
        evalCandidate(initial(q,1),initial(q,2),'initial');
    end
    for gen=1:cfg.search.proposals
        [K,alpha]=proposal(cache,gen,cfg);
        evalCandidate(K,alpha,'four_leader_proposal');
    end
    validIds=find([cache.valid]);
    if isempty(validIds)
        locked=appendRow(locked,table(f,NaN,NaN,false,string('no_feasible_train_candidate'),...
            'VariableNames',{'rollHz','K','alpha','trainFeasible','reason'}));
        fprintf(2,'%.0f Hz训练候选均无效；跳过锁定验证，保留失败日志。\n',f);
        persist();continue;
    end
    [~,loc]=min([cache(validIds).fitness]);win=cache(validIds(loc));
    locked=appendRow(locked,table(f,win.K,win.alpha,true,string('selected_on_W10_Z1_train_only'),...
        'VariableNames',{'rollHz','K','alpha','trainFeasible','reason'}));persist();
    fprintf('%.0f Hz锁定 K=%d alpha=%.4f，训练目标 %.5g\n',f,win.K,win.alpha,win.fitness);
    % 锁定后才读取验证工况。验证结果不可再反馈修改搜索参数。
    for j=1:numel(checkFiles)
        try
            dd=loadCase(checkFiles{j},f,j==2,cfg);
            chosen=solve(dd,win.K,win.alpha,true,'locked_check',cfg);
            solver=appendRow(solver,chosen.solver);
            checks=appendRow(checks,summaryRow(f,checkFiles{j},'adaptive_cond',chosen));
            % 同初始化的旧冻结K8参数作为公平对照，只有参数不同。
            if ~(win.K==8 && abs(win.alpha-154.157959862)<1e-4)
                base=solve(dd,8,154.157959862,true,'baseline_cond',cfg);
                solver=appendRow(solver,base.solver);
                checks=appendRow(checks,summaryRow(f,checkFiles{j},'frozen_cond',base));
            end
        catch ME
            fprintf(2,'检查文件 %s 不可用：%s\n',checkFiles{j},ME.message);
        end
        persist();
    end
end
writeDecision(cfg,locked,diagTbl,checks);
fprintf('E08-C结束。日志：%s\n',cfg.outputDir);

    function out=evalCandidate(K,alpha,source)
        K=round(K);alpha=round(alpha*1e4)/1e4;
        fprintf('[候选 %s] f=%.0f K=%d alpha=%.4f\n',source,f,K,alpha);
        scores=zeros(2,1);a=zeros(2,1);s=zeros(2,1);minRet=zeros(2,1);maxPhase=zeros(2,1);
        ok=true;
        for cc=1:2
            r=solve(tr{cc},K,alpha,true,['train_' source],cfg);
            solver=appendRow(solver,r.solver);persist();
            if ~r.solver.valid || ~r.metrics.sidebandPass,ok=false;end
            scores(cc)=r.fitness;a(cc)=r.metrics.nonTargetSuppressionDb(1);
            s(cc)=r.metrics.nonTargetSuppressionDb(2);
            minRet(cc)=r.metrics.minSidebandRetentionPct;
            maxPhase(cc)=r.metrics.maxPhaseShiftDeg;
        end
        if ok,fitness=mean(scores);else,fitness=Inf;end
        cache(end+1)=struct('K',K,'alpha',alpha,'fitness',fitness,'noiseA',mean(a),...
            'noiseS',mean(s),'minRetention',min(minRet),'maxPhase',max(maxPhase),'valid',ok); %#ok<AGROW>
        train=appendRow(train,table(f,K,alpha,string(source),ok,fitness,mean(a),mean(s),...
            min(minRet),max(maxPhase),...
            'VariableNames',{'rollHz','K','alpha','source','feasible','fitness',...
            'noiseProxyA_dB','noiseProxyS_dB','minRetentionPct','maxPhaseDeg'}));
        persist();out=cache(end);
    end
    function persist()
        writeIfNotEmpty(solver,fullfile(cfg.outputDir,'E08C_solver.csv'));
        writeIfNotEmpty(train,fullfile(cfg.outputDir,'E08C_train_search.csv'));
        writeIfNotEmpty(locked,fullfile(cfg.outputDir,'E08C_locked_params.csv'));
        writeIfNotEmpty(checks,fullfile(cfg.outputDir,'E08C_holdout_metrics.csv'));
        writeIfNotEmpty(diagTbl,fullfile(cfg.outputDir,'E08C_init_diagnostic.csv'));
    end
end

function r=solve(d,K,alpha,conditioned,role,cfg)
start=tic;
[modes,centers,info]=gyro_mvmd_conditioned(d.X,cfg.fs,K,alpha,...
    cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau,d.froll,conditioned);
valid=logical(info.converged && info.AC_reconRelative<=cfg.mvmd.maxAcReconstructionLoss);
fitness=Inf;
m=struct('sidebandPass',false,'nonTargetSuppressionDb',[NaN;NaN],...
    'minSidebandRetentionPct',NaN,'maxPhaseShiftDeg',NaN);
if valid
    local=cfg;local.rollHz=d.froll;
    [Y,~]=gyro_reconstruct(modes,centers,d.X,local,'paper');
    [fitness,diag]=E07D_score(d.X,Y,modes,d.t,d.frequenciesHz,d.sidebandHz,local);
    m=diag.metrics;
    valid=logical(isfinite(fitness));
end
r.fitness=fitness;r.metrics=m;
r.solver=table(string(d.fileName),d.froll,string(role),K,alpha,logical(conditioned),...
    logical(info.converged),info.iterations,info.delta,info.centerDeltaHz,...
    info.AC_reconRelative,valid,toc(start),...
    'VariableNames',{'fileName','rollHz','role','K','alpha','frequencyInitialized',...
    'converged','iterations','delta','centerDeltaHz','ACreconRelative','valid','seconds'});
end

function d=loadCase(fileName,froll,isZ,cfg)
files=[dir(fullfile(cfg.dataRoot,'**',fileName));dir(fullfile(cfg.dataRoot,fileName))];
assert(numel(files)>=1,'找不到文件：%s',fileName);
paths=unique(fullfile({files.folder},{files.name}));
assert(numel(paths)==1,'同名文件有多个路径，避免不明确的数据来源：%s',fileName);
record=read_gyro_txt(paths{1});
assert(record.invalidRows==0 && record.headerBlockRuleMatches,'文件存在异常记录或头信息：%s',fileName);
N=size(record.data,1);
iStart=1+round(cfg.trimStartSec*cfg.fs);
iStop=min(N-round(cfg.trimEndSec*cfg.fs),iStart+round(cfg.maxWindowSeconds*cfg.fs)-1);
assert(iStop>iStart+2*cfg.fs,'有效数据不足：%s',fileName);
ix=iStart:iStop;
d.fileName=fileName;d.froll=froll;
d.X=double(record.data(ix,[4 7]));d.t=(0:numel(ix)-1)'/cfg.fs;
if isZ
    fm=1;if contains(upper(fileName),'-2.'),fm=2;end
    d.frequenciesHz=[froll-fm froll froll+fm];
    d.sidebandHz=[froll-fm froll+fm];
else
    d.frequenciesHz=froll;d.sidebandHz=froll;
end
assert(all(d.frequenciesHz>0),'物理频率分量必须为正');
end

function [K,alpha]=proposal(cache,gen,cfg)
% FVIM风格四领先候选引导：不是论文作者精确优化器源码。
lo=[cfg.search.Kmin log10(cfg.search.alphaMin)];
hi=[cfg.search.Kmax log10(cfg.search.alphaMax)];
scores=[cache.fitness];[~,ids]=sort(scores);
leaders=[ [cache(ids(1:4)).K]', log10([cache(ids(1:4)).alpha]') ];
base=mean(leaders,1);a=2*(1-gen/(cfg.search.proposals+1));
for attempt=1:30
    candidates=zeros(4,2);
    for q=1:4
        A=a*(2*rand(1,2)-1);C=2*rand(1,2);
        candidates(q,:)=leaders(q,:)-A.*abs(C.*leaders(q,:)-base);
    end
    p=mean(candidates,1)+0.04*randn(1,2).*(hi-lo);
    p=max(lo,min(hi,p));
    K=round(p(1));alpha=round((10^p(2))*1e4)/1e4;
    if ~any([cache.K]==K & abs([cache.alpha]-alpha)<1e-5),return;end
end
error('E08C:duplicate','30次尝试均产生重复候选，请检查参数搜索空间');
end

function row=summaryRow(f,fileName,method,r)
row=table(f,string(fileName),string(method),r.solver.K,r.solver.alpha,...
    logical(r.solver.valid),r.fitness,r.metrics.nonTargetSuppressionDb(1),...
    r.metrics.nonTargetSuppressionDb(2),r.metrics.minSidebandRetentionPct,...
    r.metrics.maxPhaseShiftDeg,logical(r.metrics.sidebandPass),...
    'VariableNames',{'rollHz','fileName','method','K','alpha','valid','fitness',...
    'noiseProxyA_dB','noiseProxyS_dB','minRetentionPct','maxPhaseDeg','fidelityPass'});
end

function writeIfNotEmpty(T,p)
if ~isempty(T),writetable(T,p);end
end

function T=appendRow(T,row)
if isempty(T),T=row;else,T=[T;row];end
end

function writeDecision(cfg,locked,diagnostic,holdout)
path=fullfile(cfg.outputDir,'E08C_decision.txt');fid=fopen(path,'w','n','UTF-8');
assert(fid>0,'无法写入决策文件');cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'E08-C：已知滚转频率条件化初始化 + FVIM风格(K,alpha)联合参数优化可行性试验。\n');
fprintf(fid,'绝不是论文作者FVIM精确源码，不是3-21Hz全面跨频率有效性证明。\n');
fprintf(fid,'本轮仅训练Wf-10及Zf-1；Wf-20、Zf-2只在参数锁定后一次检查，不得反向调参。\n');
fprintf(fid,'W9-20已在历史研究中看过，不能称为全新盲测。\n');
for k=1:height(locked)
    fprintf(fid,'%gHz: K=%g alpha=%g 训练可行=%d\n',locked.rollHz(k),locked.K(k),locked.alpha(k),locked.trainFeasible(k));
end
for k=1:height(diagnostic)
    fprintf(fid,'初始化检查 %gHz init=%d valid=%d iter=%d delta=%g centerHz=%g AC=%g\n',...
        diagnostic.rollHz(k),diagnostic.frequencyInitialized(k),diagnostic.valid(k),...
        diagnostic.iterations(k),diagnostic.delta(k),diagnostic.centerDeltaHz(k),diagnostic.ACreconRelative(k));
end
if ~isempty(holdout)
    fprintf(fid,'完成了%d条锁定参数后的W/Z检查记录；以CSV确认收敛、保真与代理改进。\n',height(holdout));
else
    fprintf(fid,'无可用锁定后检查，需核实训练候选收敛或数据路径。\n');
end
fprintf(fid,'下一步：若3/9Hz训练和锁定后检查均有可行解，则扩展6/12/21Hz，并最后进行实际角速度和姿态RMSE验证。否则优先定位数值原因；只有合理修复仍不可行时才转单一3Hz精度优化。\n');
end
