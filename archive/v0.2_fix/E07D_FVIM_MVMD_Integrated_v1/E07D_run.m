function E07D_run()
%E07D_RUN 真实Z3-1训练：双通道连续模态重构 + 四向量FVIM风格双目标寻优。
% 主对照：原熵代理+硬重构 vs 新联合目标+连续成对重构。
% 额外四格消融：相同旧/新参数下分别硬重构与软重构，分离优化和重构效应。
% 数据边界：Z3-1已用于既往开发，训练结论不能称为独立泛化。
clc;close all;
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
cfg=E07D_config();
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
E07D_selftest(cfg);
assert(isfile(cfg.inputTrain),'缺少训练MAT文件：%s',cfg.inputTrain);
d=load(cfg.inputTrain,'X','t','baseline');
assert(isfield(d,'X')&&isfield(d,'t'),'训练MAT中缺少X或t');
X=double(d.X);t=double(d.t(:));
assert(size(X,2)==2 && size(X,1)==numel(t) && all(isfinite(X(:))));
assert(abs(median(diff(t))-1/cfg.fs)<1e-5,'训练时间轴采样率不符');
% 仅训练段评估，外部未见频率不进入参数选择和权重选择。
rng(cfg.search.randomSeed,'twister');
init=[cfg.search.initialK(:),cfg.search.initialAlpha(:)];
assert(size(init,1)==4,'应有4个共同初始候选');
cache=struct('key',{},'K',{},'alpha',{},'source',{},'valid',{}, ...
 'entropyFitness',{},'jointFitness',{},'info',{},'hard',{},'soft',{}, ...
 'hardMetrics',{},'softMetrics',{},'hardDetails',{},'softDetails',{},'modes',{},'centers',{});
trace=table();solverRows=table();
fprintf('E07-D | Z3-1训练 + 成对IMF连续重构 + 四向量FVIM式寻优\n');
fprintf('总新增分解最多%d次，单次迭代上限%d；不使用Z3-2调参。\n',...
 cfg.search.maxUniqueSolves,cfg.mvmd.maxIterations);
for i=1:4
    evaluate(init(i,1),init(i,2),'shared');
end
oldPop=init;newPop=init;
oldScore=nan(4,1);newScore=nan(4,1);
for i=1:4
    a=lookup(init(i,1),init(i,2));
    oldScore(i)=a.entropyFitness;
    newScore(i)=a.jointFitness;
end
for g=1:cfg.search.offspringPerMethod
    [oldK,oldAlpha]=proposal(oldPop,oldScore,g);
    a=evaluate(oldK,oldAlpha,'entropy');
    [~,worst]=max(oldScore);
    if a.entropyFitness<oldScore(worst)
        oldPop(worst,:)=[a.K,a.alpha];oldScore(worst)=a.entropyFitness;
    end
    [newK,newAlpha]=proposal(newPop,newScore,g);
    b=evaluate(newK,newAlpha,'joint');
    [~,worst]=max(newScore);
    if b.jointFitness<newScore(worst)
        newPop(worst,:)=[b.K,b.alpha];newScore(worst)=b.jointFitness;
    end
    [~,io]=min(oldScore);[~,in]=min(newScore);
    trace=[trace;table(g,oldPop(io,1),oldPop(io,2),oldScore(io),...
       newPop(in,1),newPop(in,2),newScore(in),numel(cache),...
       'VariableNames',{'generation','entropyK','entropyAlpha','entropyFitness',...
       'jointK','jointAlpha','jointFitness','uniqueSolves'})]; %#ok<AGROW>
    saveProgress();
    fprintf('gen %02d/%02d | unique solves=%d | old %.6g | joint %.6g\n',...
       g,cfg.search.offspringPerMethod,numel(cache),oldScore(io),newScore(in));
end
[~,io]=min(oldScore);[~,in]=min(newScore);
old=lookup(oldPop(io,1),oldPop(io,2));
new=lookup(newPop(in,1),newPop(in,2));
assert(old.valid && new.valid,'未找到可用的收敛候选，检查求解器输出。');
% 同一共同候选池内的事后全池排名，仅诊断；不会冒称各自算法曾访问全部候选。
[~,ioAll]=min([cache.entropyFitness]);[~,inAll]=min([cache.jointFitness]);
pool=table(string({'entropy_globalpool';'joint_globalpool'}),...
    [cache(ioAll).K;cache(inAll).K],[cache(ioAll).alpha;cache(inAll).alpha],...
    [cache(ioAll).entropyFitness;cache(inAll).jointFitness],...
    'VariableNames',{'role','K','alpha','fitness'});
writetable(pool,fullfile(cfg.outputDir,'E07D_shared_pool_diagnostic.csv'));
% 两参数×两重构：只比较同一个缓存的MVMD分解，提供正交消融。
rows=table();
items={old,new};names={'entropy_selected','joint_selected'};
for j=1:2
    c=items{j};
    rows=[rows;resultRow('Z3-1',names{j},'hard',c.K,c.alpha,c.hardMetrics)]; %#ok<AGROW>
    rows=[rows;resultRow('Z3-1',names{j},'soft',c.K,c.alpha,c.softMetrics)]; %#ok<AGROW>
end
writetable(rows,fullfile(cfg.outputDir,'E07D_train_ablation.csv'));
writetable(old.hardDetails,fullfile(cfg.outputDir,'E07D_old_IMF.csv'));
writetable(new.softDetails,fullfile(cfg.outputDir,'E07D_new_IMF.csv'));
% 两参数若一致是有价值的诊断，但不可虚构优化贡献。
params=table(string({'entropy';'joint'}),[old.K;new.K],[old.alpha;new.alpha],...
  [old.entropyFitness;new.jointFitness],...
  'VariableNames',{'objective','K','alpha','fitness'});
writetable(params,fullfile(cfg.outputDir,'E07D_locked_params.csv'));
save(fullfile(cfg.outputDir,'E07D_train_outputs.mat'),...
 'X','t','params','rows','trace','cfg','-v7');
try
    plotTraining(X,old.hard,new.soft,cfg);
catch err
    fprintf(2,'绘图略过（不影响计算）：%s\n',err.message);
end
% 检查记录：新策略是否同时保持动态边带、并降低两个通道残差代理。
passSide=new.softMetrics.sidebandPass;
gainA=new.softMetrics.nonTargetSuppressionDb(1)-old.hardMetrics.nonTargetSuppressionDb(1);
gainS=new.softMetrics.nonTargetSuppressionDb(2)-old.hardMetrics.nonTargetSuppressionDb(2);
writeDecision(old,new,passSide,gainA,gainS,cfg);
% 锁定权重和K-alpha后，若提供完整历史TXT目录，才运行真正跨频率的留出评价。
if ~isempty(cfg.dataRoot)
    E07D_holdout(old,new,cfg);
else
    fprintf('未设置dataRoot：已输出训练结果，未运行外部留出数据。\n');
end
fprintf('E07-D结束，输出：%s\n',cfg.outputDir);

    function item=evaluate(K,alpha,role)
        K=round(K);alpha=max(cfg.search.alphaMin,min(cfg.search.alphaMax,alpha));
        alpha=round(alpha*1e4)/1e4;
        key=sprintf('K%d_A%.4f',K,alpha);
        ix=find(strcmp({cache.key},key),1);
        if ~isempty(ix),item=cache(ix);return;end
        assert(numel(cache)<cfg.search.maxUniqueSolves,'已达MVMD分解预算上限');
        fprintf('[%02d/%02d] %s: K=%d, alpha=%.4f\n',...
          numel(cache)+1,cfg.search.maxUniqueSolves,role,K,alpha);
        timer=tic;
        [modes,centers,info]=gyro_mvmd(X,cfg.fs,K,alpha,...
            cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
        elapsed=toc(timer);
        valid=logical(info.converged && ...
            info.AC_reconRelative<=cfg.joint.maxAcReconstructionLoss);
        scoreOld=Inf;scoreNew=Inf;hard=[];soft=[];hm=[];sm=[];
        hd=table();sd=table();
        if valid
            [hard,hd]=gyro_reconstruct(modes,centers,X,cfg,'paper');
            [soft,sd]=E07D_soft_reconstruct(modes,centers,X,cfg);
            scoreOld=E07D_entropy(modes);
            [scoreNew,diagnostic]=E07D_score(X,soft,modes,t,...
                cfg.trainFrequenciesHz,cfg.trainSidebandsHz,cfg);
            sm=diagnostic.metrics;
            hm=E07B_candidate_metrics(X,hard,t,...
                cfg.trainFrequenciesHz,cfg.trainSidebandsHz,cfg);
            valid=isfinite(scoreOld)&&isfinite(scoreNew);
        end
        if ~valid,scoreOld=Inf;scoreNew=Inf;end
        item=struct('key',key,'K',K,'alpha',alpha,'source',role,...
          'valid',valid,'entropyFitness',scoreOld,'jointFitness',scoreNew,...
          'info',info,'hard',hard,'soft',soft,'hardMetrics',hm,'softMetrics',sm,...
          'hardDetails',hd,'softDetails',sd,'modes',modes,'centers',centers);
        cache(end+1)=item;
        % 每个候选耗时写入轻量CSV，遇中断可以判断停止位置。
        rec=table(string(key),string(role),K,alpha,valid,logical(info.converged),...
           info.iterations,info.AC_reconRelative,scoreOld,scoreNew,elapsed,...
           'VariableNames',{'key','role','K','alpha','valid','converged',...
           'iterations','ACreconRelative','entropyFitness','jointFitness','seconds'});
        p=fullfile(cfg.outputDir,'E07D_solvers.csv');
        solverRows=[solverRows;rec]; %#ok<AGROW>
        writetable(solverRows,p);
    end
    function a=lookup(K,alpha)
        key=sprintf('K%d_A%.4f',round(K),round(alpha*1e4)/1e4);
        ix=find(strcmp({cache.key},key),1);
        assert(~isempty(ix),'候选缓存丢失');a=cache(ix);
    end
    function [K,alpha]=proposal(pop,score,g)
        lo=[cfg.search.Kmin log10(cfg.search.alphaMin)];
        hi=[cfg.search.Kmax log10(cfg.search.alphaMax)];
        u=[pop(:,1) log10(pop(:,2))];
        [~,ids]=sort(score);leaders=u(ids(1:4),:);
        a=2*(1-g/(cfg.search.offspringPerMethod+1));
        base=mean(u,1);trial=base;
        for q=1:20
            moves=zeros(4,2);
            for j=1:4
                A=a*(2*rand(1,2)-1);C=2*rand(1,2);
                dist=abs(C.*leaders(j,:)-base);
                moves(j,:)=leaders(j,:)-A.*dist;
            end
            trial=mean(moves,1);
            if mod(g,3)==0
                trial=trial+0.09*randn(1,2).*(hi-lo);
            end
            trial=max(lo,min(hi,trial));
            K=round(trial(1));alpha=round(10^trial(2)*1e4)/1e4;
            key=sprintf('K%d_A%.4f',K,alpha);
            if ~any(strcmp({cache.key},key)),return;end
            % 重复候选时确定性随机探索，避免多次求相同分解。
            base=lo+rand(1,2).*(hi-lo);
        end
        error('E07D:duplicateCandidates','连续20次产生重复候选；停止而不扩大计算预算');
    end
    function saveProgress()
        writetable(trace,fullfile(cfg.outputDir,'E07D_search_trace.csv'));
    end
end

function row=resultRow(caseName,selectedBy,reconstruction,K,alpha,m)
row=table(string(caseName),string(selectedBy),string(reconstruction),K,alpha,...
  m.nonTargetSuppressionDb(1),m.nonTargetSuppressionDb(2),...
  m.minSidebandRetentionPct,m.maxSidebandRetentionPct,...
  m.maxPhaseShiftDeg,m.maxRelativePhaseShiftDeg,logical(m.sidebandPass),...
  'VariableNames',{'caseName','selectedBy','reconstruction','K','alpha',...
  'noiseProxyDbA','noiseProxyDbS','minRetentionPct','maxRetentionPct',...
  'maxPhaseShiftDeg','maxRelativeShiftDeg','sidebandPass'});
end

function plotTraining(X,oldY,newY,cfg)
N=size(X,1);f=(0:floor(N/2))'*cfg.fs/N;
fig=figure('Visible','off','Color','w','Position',[100 100 980 680]);
for c=1:2
  subplot(2,1,c);hold on;
  spectrum=@(x) 2*abs(fft(x-mean(x)))/numel(x);
  x=spectrum(X(:,c));o=spectrum(oldY(:,c));n=spectrum(newY(:,c));
  plot(f,x(1:numel(f)),'Color',[.6 .6 .6]);
  plot(f,o(1:numel(f)),'Color',[.85 .4 .1]);
  plot(f,n(1:numel(f)),'Color',[.1 .45 .8]);
  grid on;xlim([0 15]);ylabel(sprintf('Channel %d amplitude',c));
  if c==1,title('E07-D Z3-1 TRAIN: old-hard vs new-soft');
    legend('Raw','Entropy + hard','Joint + soft');end
end
xlabel('Frequency (Hz)');
saveas(fig,fullfile(cfg.outputDir,'E07D_train_spectrum.png'));
close(fig);
end

function writeDecision(old,new,side,gainA,gainS,cfg)
fid=fopen(fullfile(cfg.outputDir,'E07D_decision.txt'),'w','n','UTF-8');
assert(fid>0,'无法写入E07D_decision');cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'E07-D 真实Z3-1 TRAIN；FVIM风格四引导向量工程复现，非论文作者精确源码。\n');
fprintf(fid,'旧目标K=%d alpha=%.5f; 新目标K=%d alpha=%.5f。\n',...
 old.K,old.alpha,new.K,new.alpha);
fprintf(fid,'相对旧目标+硬重构：新目标+软重构额外代理抑噪 A=%.5f,S=%.5f dB。\n',gainA,gainS);
fprintf(fid,'训练边带保真是否达标：%d（99..101%%，相位<0.2度）。\n',side);
if side && gainA>0 && gainS>0
    fprintf(fid,'训练可行性支持继续外部工况检查；不是独立验证通过。\n');
else
    fprintf(fid,'未同时满足双路净改善与保真，仍检查消融以定位优化/重构贡献；不改留出数据上的参数。\n');
end
fprintf(fid,'固定重构规则：csWidth=%g, sampWidth=%g, mixExponent=%g。\n',...
 cfg.soft.csWidth,cfg.soft.sampWidth,cfg.soft.mixExponent);
fprintf(fid,'仅为单频拟合残差代理，不是真值SNR、角速度RMSE或姿态角RMSE。\n');
end
