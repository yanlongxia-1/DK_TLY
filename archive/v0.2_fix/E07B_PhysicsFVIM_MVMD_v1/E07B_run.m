function E07B_run()
% E07-B｜低预算、可缓存、双目标四优胜者引导FVIM风格搜索。
% 与论文作者的FVIM源码不等同。只在训练Z3-1上选参数。
clc; close all;
cfg=E07B_config();
assert(isfile(cfg.inputTrain),'找不到训练MAT数据: %s',cfg.inputTrain);
assert(isfile(cfg.inputCheck),'找不到检查MAT数据: %s',cfg.inputCheck);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
E07B_selftest();
rng(cfg.search.randomSeed,'twister');
tr=loadInput(cfg.inputTrain);
assert(numel(tr.t)==size(tr.X,1),'训练时间与样本长度不匹配');
% 此次搜索只使用 Z3-1 的原始双通道测量及先前冻结MVMD输出。
baseMetrics=E07B_candidate_metrics(tr.X,tr.baseline,tr.t, ...
    cfg.trainComponentsHz,cfg.trainSidebandHz,cfg);
cache=struct('key',{},'K',{},'alpha',{},'valid',{},'oldFit',{},'jointFit',{}, ...
    'metrics',{},'Y',{},'info',{},'modes',{},'centers',{},'details',{});
solverRows=table();candidateRows=table();searchTrace=table();
maxN=cfg.search.maxUniqueTrainSolves;
init=[cfg.search.initialK(:),cfg.search.initialAlpha(:)];
assert(size(init,1)==4,'四优胜者初始群体恰好4个');

fprintf('\nE07-B | dual-channel FVIM-style joint-objective pilot\n');
fprintf('训练: Z3-1 | 四个共同初始点 | 每目标四次有限更新 | 最多%d次不同MVMD求解\n',maxN);
for k=1:4
    evalPoint(init(k,1),init(k,2),'shared',0); %#ok<NASGU>
end
oldPop=init;newPop=init;
oldScore=zeros(4,1);newScore=zeros(4,1);
for i=1:4
    z=getByParam(init(i,1),init(i,2));
    oldScore(i)=z.oldFit;
    newScore(i)=z.jointFit;
end
for gen=1:cfg.search.offSpringPerMethod
    % 同一随机流顺序、相同四优胜者结构，但按各自目标排序并产生子代。
    [k1,a1]=suggest(oldPop,oldScore,gen,cfg);
    z1=evalPoint(k1,a1,'entropy',gen);
    [worstScore,worstIdx]=max(oldScore);
    if z1.oldFit<worstScore
        oldPop(worstIdx,:)=[k1 a1];oldScore(worstIdx)=z1.oldFit;
    end
    [k2,a2]=suggest(newPop,newScore,gen,cfg);
    z2=evalPoint(k2,a2,'joint',gen);
    [worstScore,worstIdx]=max(newScore);
    if z2.jointFit<worstScore
        newPop(worstIdx,:)=[k2 a2];newScore(worstIdx)=z2.jointFit;
    end
    [bestOld,oldIdx]=min(oldScore);[bestJoint,jointIdx]=min(newScore);
    searchTrace=[searchTrace;table(gen,oldPop(oldIdx,1),oldPop(oldIdx,2),bestOld,...
        newPop(jointIdx,1),newPop(jointIdx,2),bestJoint, ...
        'VariableNames',{'generation','entropyBestK','entropyBestAlpha','entropyBestFitness',...
        'jointBestK','jointBestAlpha','jointBestFitness'})]; %#ok<AGROW>
    fprintf('  gen %d/%d: entropy K=%d a=%.5g | joint K=%d a=%.5g | cached solves=%d\n',...
        gen,cfg.search.offSpringPerMethod,oldPop(oldIdx,1),oldPop(oldIdx,2),...
        newPop(jointIdx,1),newPop(jointIdx,2),numel(cache));
end
assert(numel(cache)<=maxN,'超出训练MVMD最大求解预算');
[~,io]=min(oldScore); [~,ij]=min(newScore);
old=getByParam(oldPop(io,1),oldPop(io,2));
joint=getByParam(newPop(ij,1),newPop(ij,2));
allTrain=table();
allTrain=[allTrain; metricRow('Z3-1','frozen_MVMD',8,cfg.mvmd.alphaBaseline,baseMetrics,baseMetrics,cfg)];
if isfinite(old.oldFit)
    allTrain=[allTrain;metricRow('Z3-1','entropy_engineering',old.K,old.alpha,old.metrics,baseMetrics,cfg)];
end
if isfinite(joint.jointFit)
    allTrain=[allTrain;metricRow('Z3-1','joint_physics',joint.K,joint.alpha,joint.metrics,baseMetrics,cfg)];
end
saveTables();
trainPass=false;
if isfinite(joint.jointFit)
    newExtra=joint.metrics.nonTargetSuppressionDb-baseMetrics.nonTargetSuppressionDb;
    againstOld=false;
    if isfinite(old.oldFit)
        oldExtra=old.metrics.nonTargetSuppressionDb-baseMetrics.nonTargetSuppressionDb;
        againstOld=all(newExtra-oldExtra>=cfg.pass.minAgainstEntropyDb) || ~old.metrics.sidebandPass;
    else
        againstOld=true;
    end
    trainPass=joint.metrics.sidebandPass && ...
        all(newExtra>=cfg.pass.minExtraSuppressionDb) && againstOld;
end
if ~trainPass
    decision=['STOP: joint FVIM-style objective did not meet preset train conditions. ',...
        'Do NOT inspect Z3-2 to retune. Record failure and move on.'];
    writeDecision(decision);
    fprintf('\n%s\n',decision);
    return;
end
% 满足训练约束，固定结果后才使用Z3-2。本数据在E06A/B阶段已经被看过，非项目全新盲测。
te=loadInput(cfg.inputCheck);
baseCheck=E07B_candidate_metrics(te.X,te.baseline,te.t, ...
    cfg.checkComponentsHz,cfg.checkSidebandHz,cfg);
checkRows=metricRow('Z3-2','frozen_MVMD',8,cfg.mvmd.alphaBaseline,baseCheck,baseCheck,cfg);
checkParams=[old.K old.alpha; joint.K joint.alpha];
checkMethod={'entropy_engineering','joint_physics'};
checkCache=struct('key',{},'Y',{},'valid',{});
for i=1:2
    K=checkParams(i,1);a=checkParams(i,2);key=keyOf(K,a);
    idx=find(strcmp({checkCache.key},key),1);
    if isempty(idx)
        fprintf('锁定参数检查 Z3-2：%s K=%d a=%.6g\n',checkMethod{i},K,a);
        tic;[modes,centerHz,info]=gyro_mvmd(te.X,cfg.fs,K,a, ...
            cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
        elapsed=toc;
        valid=logical(info.converged && info.AC_reconRelative<=cfg.joint.maxAcReconstructionLoss);
        if valid
            [Y,~]=gyro_reconstruct(modes,centerHz,te.X,cfg,'paper');
        else
            Y=nan(size(te.X));
        end
        checkCache(end+1)=struct('key',key,'Y',Y,'valid',valid); %#ok<AGROW>
        solverRows=[solverRows;table(string('Z3-2'),string('check'),K,a,...
            logical(info.converged),info.iterations,info.AC_reconRelative,elapsed,...
            'VariableNames',{'caseName','role','K','alpha','converged','iterations',...
            'acReconRelative','elapsedSeconds'})]; %#ok<AGROW>
        idx=numel(checkCache);
    end
    if checkCache(idx).valid
        m=E07B_candidate_metrics(te.X,checkCache(idx).Y,te.t, ...
            cfg.checkComponentsHz,cfg.checkSidebandHz,cfg);
        checkRows=[checkRows;metricRow('Z3-2',checkMethod{i},K,a,m,baseCheck,cfg)]; %#ok<AGROW>
    else
        fprintf(2,'锁定参数在Z3-2未收敛或重构误差超过限制：%s\n',checkMethod{i});
    end
end
writetable(checkRows,fullfile(cfg.outputDir,'E07B_locked_check.csv'));
writetable(solverRows,fullfile(cfg.outputDir,'E07B_solver.csv'));
checkJoint=checkRows(strcmp(checkRows.method,'joint_physics'),:);
checkEntropy=checkRows(strcmp(checkRows.method,'entropy_engineering'),:);
if height(checkJoint)==1 && height(checkEntropy)==1 && ...
   checkJoint.sidebandPass && ...
   checkJoint.extraNoiseDbA>=cfg.pass.minExtraSuppressionDb && ...
   checkJoint.extraNoiseDbS>=cfg.pass.minExtraSuppressionDb && ...
   min(checkJoint.extraNoiseDbA-checkEntropy.extraNoiseDbA,...
       checkJoint.extraNoiseDbS-checkEntropy.extraNoiseDbS)>=cfg.pass.minAgainstEntropyDb
    decision='PROVISIONAL: joint objective passed locked Z3-2 check. NOT independent blind proof; reserve NEW unseen conditions.';
else
    decision='STOP after locked Z3-2 check: no reproducible gain; do not tune on Z3-2.';
end
writeDecision(decision);
fprintf('\n%s\n',decision);

    function z=evalPoint(K,alpha,role,generation)
        K=round(K);alpha=double(alpha);
        key=keyOf(K,alpha);
        ii=find(strcmp({cache.key},key),1);
        if ~isempty(ii)
            z=cache(ii);
            candidateRows=[candidateRows; candidateRow(z,role,generation,1)];
            return;
        end
        assert(numel(cache)<maxN,'参数候选分解预算耗尽');
        fprintf('[TRAIN %02d/%02d] %s gen=%d K=%d alpha=%.7g\n', ...
            numel(cache)+1,maxN,role,generation,K,alpha);
        tic;
        [modes,centers,info]=gyro_mvmd(tr.X,cfg.fs,K,alpha, ...
            cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
        elapsed=toc;
        valid=logical(info.converged && info.AC_reconRelative<=cfg.joint.maxAcReconstructionLoss);
        oldFit=Inf;jointFit=Inf;m=[];Y=[];details=[];
        if valid
            [Y,details]=gyro_reconstruct(modes,centers,tr.X,cfg,'paper');
            [oldFit,jointFit,~,m]=E07B_objective(tr.X,Y,modes,tr.t, ...
                cfg.trainComponentsHz,cfg.trainSidebandHz,cfg);
            valid=isfinite(oldFit)&&isfinite(jointFit);
            if ~valid,oldFit=Inf;jointFit=Inf;end
        end
        z=struct('key',key,'K',K,'alpha',alpha,'valid',valid,'oldFit',oldFit, ...
            'jointFit',jointFit,'metrics',m,'Y',Y,'info',info, ...
            'modes',[],'centers',centers,'details',details);
        cache(end+1)=z;
        solverRows=[solverRows;table(string('Z3-1'),string(role),K,alpha,...
            logical(info.converged),info.iterations,info.AC_reconRelative,elapsed,...
            'VariableNames',{'caseName','role','K','alpha','converged',...
            'iterations','acReconRelative','elapsedSeconds'})]; %#ok<AGROW>
        candidateRows=[candidateRows; candidateRow(z,role,generation,0)];
    end
    function z=getByParam(K,alpha)
        ii=find(strcmp({cache.key},keyOf(K,alpha)),1);
        assert(~isempty(ii),'尚未计算请求的候选');z=cache(ii);
    end
    function row=candidateRow(z,role,generation,cached)
        if z.valid
            ab=z.metrics.maxAmplitudeDeviationPct;
            pd=z.metrics.maxPhaseShiftDeg;
            rp=z.metrics.maxRelativePhaseShiftDeg;
            nd=z.metrics.nonTargetSuppressionDb(:);
            sb=z.metrics.sidebandPass;
        else
            ab=NaN;pd=NaN;rp=NaN;nd=[NaN;NaN];sb=false;
        end
        row=table(string('Z3-1'),string(role),generation,logical(cached),...
            z.K,z.alpha,logical(z.valid),z.oldFit,z.jointFit,nd(1),nd(2),...
            ab,pd,rp,logical(sb),'VariableNames',{'caseName','role',...
            'generation','cached','K','alpha','valid','entropyFitness','jointFitness',...
            'nonTargetSuppressionDbA','nonTargetSuppressionDbS','maxAmplitudeDeviationPct',...
            'maxPhaseAbsDeg','maxRelativePhaseDeg','sidebandPass'});
    end
    function saveTables()
        writetable(solverRows,fullfile(cfg.outputDir,'E07B_solver.csv'));
        writetable(candidateRows,fullfile(cfg.outputDir,'E07B_train_candidates.csv'));
        writetable(searchTrace,fullfile(cfg.outputDir,'E07B_search_trace.csv'));
        writetable(allTrain,fullfile(cfg.outputDir,'E07B_train_comparison.csv'));
        % 不把MAT写入大体积的分解模态，仅保存已评价信息和最佳两组输出。
        if isfinite(old.oldFit) && isfinite(joint.jointFit)
            trainRaw=tr.X;trainFrozen=tr.baseline;
            trainEntropy=old.Y;trainJoint=joint.Y;trainTime=tr.t;
            save(fullfile(cfg.outputDir,'E07B_train_outputs.mat'),...
                'trainRaw','trainFrozen','trainEntropy','trainJoint','trainTime','-v7');
            if istable(old.details)
                writetable(old.details,fullfile(cfg.outputDir,'E07B_entropy_IMF_modes.csv'));
            end
            if istable(joint.details)
                writetable(joint.details,fullfile(cfg.outputDir,'E07B_joint_IMF_modes.csv'));
            end
            % 频谱图为说明性材料，不参与搜索、排名或停止条件。
            try
                fig=figure('Visible','off','Color','w','Position',[100 100 960 700]);
                for c=1:2
                    subplot(2,1,c);
                    f=(0:floor(size(tr.X,1)/2))'*cfg.fs/size(tr.X,1);
                    r=localSpectrum(tr.X(:,c));
                    b=localSpectrum(tr.baseline(:,c));
                    o=localSpectrum(old.Y(:,c));
                    n=localSpectrum(joint.Y(:,c));
                    plot(f,r,'Color',[.6 .6 .6]);hold on;
                    plot(f,b,'b');plot(f,o,'Color',[.85 .45 .1]);
                    plot(f,n,'Color',[.15 .55 .3]);
                    xlim([0 15]);grid on;
                    ylabel(sprintf('Channel %d amplitude',c));
                    if c==1
                        title('Z3-1 train spectra | old entropy vs joint objective');
                        legend('raw','frozen MVMD','entropy fitness','joint physics fitness');
                    end
                end
                xlabel('Frequency (Hz)');
                saveas(fig,fullfile(cfg.outputDir,'E07B_train_spectrum.png'));
                close(fig);
            catch plotErr
                fprintf(2,'绘图失败（不影响表格和实验判定）：%s\n',plotErr.message);
            end
        end
    end
    function writeDecision(message)
        fid=fopen(fullfile(cfg.outputDir,'E07B_decision.txt'),'w','n','UTF-8');
        assert(fid>0,'无法写入decision文件');
        fprintf(fid,'%s\n',message);
        fprintf(fid,'Old best (engineering entropy-only): K=%d alpha=%.9g fitness=%.9g\n',...
            old.K,old.alpha,old.oldFit);
        fprintf(fid,'New best (dual-channel physics objective): K=%d alpha=%.9g fitness=%.9g\n',...
            joint.K,joint.alpha,joint.jointFit);
        fprintf(fid,'Unique train MVMD solves: %d, max: %d. Locked Z3-2 is NOT unseen blind test.\n',numel(cache),maxN);
        fprintf(fid,'No ground-truth signal SNR, body Ωx/Ωy RMSE, or attitude RMSE demonstrated.\n');
        fclose(fid);
    end
end

function key=keyOf(K,alpha)
key=sprintf('K%d_a%.5f',round(K),alpha);
end
function d=loadInput(fileName)
d=load(fileName,'X','baseline','t');
assert(all(isfield(d,{'X','baseline','t'})),'MAT缺少 X / baseline / t');
assert(isequal(size(d.X),size(d.baseline)) && size(d.X,2)==2,'MAT双通道维度无效');
assert(all(isfinite(d.X(:))) && all(isfinite(d.baseline(:))),'输入存在NaN/Inf');
d.t=d.t(:);
assert(size(d.X,1)==numel(d.t),'时间与数据长度不同');
end
function [K,alpha]=suggest(pop,scores,gen,cfg)
% 四优胜者引导的 FVIM-style 子代提议；不是作者FVIM源码逐位复现。
[~,idx]=sort(scores,'ascend');leaders=pop(idx,:);
base=pop(mod(gen-1,4)+1,:);
lo=[cfg.search.Kmin log10(cfg.search.alphaMin)];
hi=[cfg.search.Kmax log10(cfg.search.alphaMax)];
p=[base(1),log10(base(2))];
L=[leaders(:,1),log10(leaders(:,2))];
a=1.2*(1-gen/(cfg.search.offSpringPerMethod+1));
moves=zeros(4,2);
for j=1:4
    A=a*(2*rand(1,2)-1);C=2*rand(1,2);
    moves(j,:)=L(j,:)-A.*abs(C.*L(j,:)-p);
end
p2=mean(moves,1);
% 确定性约束与量化，防止为数值上几乎相同的alpha重复求解。
p2=min(hi,max(lo,p2));
K=max(cfg.search.Kmin,min(cfg.search.Kmax,round(p2(1))));
alpha=round((10^p2(2))*1e4)/1e4; % 统一量化至1e-4，避免不同版本round(...,'significant')兼容性
alpha=min(cfg.search.alphaMax,max(cfg.search.alphaMin,alpha));
end
function row=metricRow(caseName,method,K,alpha,m,base,cfg)
extra=m.nonTargetSuppressionDb-base.nonTargetSuppressionDb;
pass=logical(m.sidebandPass && all(extra>=cfg.pass.minExtraSuppressionDb));
row=table(string(caseName),string(method),K,alpha,...
    m.nonTargetSuppressionDb(1),m.nonTargetSuppressionDb(2),...
    extra(1),extra(2),m.minSidebandRetentionPct,m.maxSidebandRetentionPct,...
    m.maxPhaseShiftDeg,m.maxRelativePhaseShiftDeg,...
    logical(m.sidebandPass),pass,...
    'VariableNames',{'caseName','method','K','alpha',...
    'nonTargetSuppressionDbA','nonTargetSuppressionDbS','extraNoiseDbA','extraNoiseDbS',...
    'minSidebandRetentionPct','maxSidebandRetentionPct','maxSidebandPhaseShiftDeg',...
    'maxRelativePhaseShiftDeg','sidebandPass','passedVsFrozen'});
end

function amp=localSpectrum(x)
N=numel(x);X=fft(x-mean(x));
amp=2*abs(X(1:floor(N/2)+1))/N;
amp(1)=amp(1)/2;
if rem(N,2)==0,amp(end)=amp(end)/2;end
end
