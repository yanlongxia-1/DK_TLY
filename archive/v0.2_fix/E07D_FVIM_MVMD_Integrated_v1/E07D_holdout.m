function E07D_holdout(old,new,cfg)
%E07D_HOLDOUT 参数锁定后跨滚转频率的真实数据检验；只用于报告，不回馈训练。
% 未上传整套TXT时，不臆造验证结果。保证用完全相同预处理和采样率。
assert(isfolder(cfg.dataRoot),'设置的原始TXT目录不存在');
rows=table();solver=table();
for i=1:numel(cfg.holdout)
    fileName=cfg.holdout{i};
    allMatches=dir(fullfile(cfg.dataRoot,'**',fileName));
    if isempty(allMatches)
        fprintf(2,'留出文件尚未找到，跳过并记录：%s\n',fileName);
        continue;
    end
    if numel(allMatches)>1
        error('E07D:ambiguousFile','留出文件%s有多个路径，必须消歧。',fileName);
    end
    filePath=fullfile(allMatches(1).folder,allMatches(1).name);
    record=read_gyro_txt(filePath);
    data=record.data;
    N=size(data,1);
    first=1+round(cfg.holdoutTrimStart*cfg.fs);
    last=N-round(cfg.holdoutTrimEnd*cfg.fs);
    assert(last-first>=2*cfg.fs,'留出文件有效数据不足2秒：%s',fileName);
    data=data(first:last,:);
    X=double(data(:,[4,7]));
    t=(0:size(X,1)-1)'/cfg.fs;
    rollHz=cfg.holdoutRollHz(i);
    modHz=cfg.holdoutMotionHz(i);
    localCfg=cfg;localCfg.rollHz=rollHz;
    if modHz>0
        frequencies=unique(sort([rollHz-modHz,rollHz,rollHz+modHz]));
        sidebands=[rollHz-modHz,rollHz+modHz];
    else
        frequencies=rollHz;
        sidebands=rollHz;
    end
    params=[old.K old.alpha;new.K new.alpha];
    selectors={'entropy','joint'};
    cacheKeys={};cacheResults={};
    for j=1:2
        K=params(j,1);alpha=params(j,2);
        key=sprintf('K%d_A%.4f',K,alpha);
        idx=find(strcmp(cacheKeys,key),1);
        if isempty(idx)
            fprintf('留出检查%s: K=%d alpha=%g\n',fileName,K,alpha);
            tic;
            [modes,centerHz,info]=gyro_mvmd(X,cfg.fs,K,alpha,...
                cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
            elapsed=toc;
            valid=logical(info.converged && ...
                 info.AC_reconRelative<=cfg.joint.maxAcReconstructionLoss);
            if valid
                [Yhard,detailsH]=gyro_reconstruct(modes,centerHz,X,localCfg,'paper');
                [Ysoft,detailsS]=E07D_soft_reconstruct(modes,centerHz,X,localCfg);
            else
                Yhard=[];Ysoft=[];detailsH=table();detailsS=table();
            end
            cacheKeys{end+1}=key; %#ok<AGROW>
            cacheResults{end+1}=struct('valid',valid,'hard',Yhard,'soft',Ysoft,...
                'hardDetails',detailsH,'softDetails',detailsS); %#ok<AGROW>
            idx=numel(cacheKeys);
            solver=[solver;table(string(fileName),K,alpha,...
               logical(info.converged),info.iterations,info.AC_reconRelative,elapsed,...
               'VariableNames',{'fileName','K','alpha','converged','iterations',...
               'ACreconRelative','elapsedSeconds'})]; %#ok<AGROW>
        end
        item=cacheResults{idx};
        if ~item.valid
            fprintf(2,'警告：留出工况%s分解不满足收敛或重构约束\n',fileName);
            continue;
        end
        for p=1:2
            if p==1,Y=item.hard;policy='hard';else,Y=item.soft;policy='soft';end
            m=E07B_candidate_metrics(X,Y,t,frequencies,sidebands,localCfg);
            rows=[rows;table(string(fileName),rollHz,string(selectors{j}),...
                string(policy),K,alpha,m.nonTargetSuppressionDb(1),...
                m.nonTargetSuppressionDb(2),m.minSidebandRetentionPct,...
                m.maxSidebandRetentionPct,m.maxPhaseShiftDeg,...
                m.maxRelativePhaseShiftDeg,logical(m.sidebandPass),...
                'VariableNames',{'fileName','rollHz','selectedBy','reconstruction',...
                'K','alpha','noiseProxyDbA','noiseProxyDbS','minRetentionPct',...
                'maxRetentionPct','maxPhaseShiftDeg','maxRelativeShiftDeg',...
                'sidebandPass'})]; %#ok<AGROW>
        end
    end
end
if ~isempty(rows)
    writetable(rows,fullfile(cfg.outputDir,'E07D_holdout_comparison.csv'));
end
if ~isempty(solver)
    writetable(solver,fullfile(cfg.outputDir,'E07D_holdout_solver.csv'));
end
end
