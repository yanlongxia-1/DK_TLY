function E08F_run()
%E08F_RUN: diagnose E08E's all-1200 nonconvergence and short Z9-3 segment.
% Does not change alpha/K, objective, tau, tolerance, calibration, or reconstruction.
% W9-60 and Z9-3 have now been inspected and are DEVELOPMENT diagnostics, not blind tests.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
cfg=E08F_config();
assert(isfolder(cfg.dataRoot),'E08F:Path','Set cfg.dataRoot in E08F_config.m');
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
windows=table();results=table();history=table();issues=table();
for fi=1:numel(cfg.filenames)
    name=cfg.filenames{fi};
    try
        path=locateOne(cfg.dataRoot,name);
        r=read_gyro_txt(path);
        assert(r.invalidRows==0 && r.headerBlockRuleMatches,'TXT invalid rows or header mismatch');
        N=size(r.data,1);
        % Previous E08E selection, recorded exactly for comparison.
        oldStart=max(1,round(cfg.testFraction(1)*N));
        oldEnd=min([N, round(cfg.testFraction(2)*N),...
             oldStart+round(cfg.testMaxSec*cfg.fs)-1]);
        % Consistent tail-window diagnostic. Preserve 0.5s startup exclusion and 5%% end gap.
        tailEnd=min(N,max(1,floor(cfg.testFraction(2)*N)));
        tailStart=max(1+round(cfg.skipStartSec*cfg.fs),...
             tailEnd-round(cfg.testMaxSec*cfg.fs)+1);
        starts=[oldStart tailStart];ends=[oldEnd tailEnd];labels={'legacy_60to95','tail_last3s'};
        nrows=max(0,ends-starts+1);
        w=table(string(name),N,starts(1),ends(1),nrows(1),starts(2),ends(2),nrows(2),...
            logical(nrows(1)>=cfg.minTestSec*cfg.fs),...
            logical(nrows(2)>=cfg.minTestSec*cfg.fs),...
            'VariableNames',{'fileName','totalRows','legacyStart','legacyEnd','legacyRows',...
            'tailStart','tailEnd','tailRows','legacyLongEnough','tailLongEnough'});
        windows=appendRows(windows,w);
        fprintf('\n%s: total=%d, legacy=%d, tail=%d samples\n',name,N,nrows(1),nrows(2));
        for wi=1:numel(labels)
            if nrows(wi)<cfg.minTestSec*cfg.fs
                issues=appendRows(issues,newIssue(name,labels{wi},'window_below_2s'));
                fprintf('  SKIP %s: %.3fs < %.3fs\n',labels{wi},nrows(wi)/cfg.fs,cfg.minTestSec);
                continue;
            end
            X=double(r.data(starts(wi):ends(wi),[4 7]));
            for k=1:size(cfg.candidates,1)
                K=cfg.candidates(k,1);alpha=cfg.candidates(k,2);
                fprintf('  %s %s K=%d alpha=%.4f ...\n',labels{wi},cfg.methods{k},K,alpha);
                ticRun=tic;
                try
                    [modes,cent,info]=gyro_mvmd_conditioned(X,cfg.fs,K,alpha,...
                        cfg.maxIterDiagnostic,cfg.tol,cfg.tau,cfg.rollHz,true);
                    secs=toc(ticRun);
                    n=info.iterations;h=info.deltaHistory;g=info.centerHistoryHz;
                    at1200=NaN;centerAt1200=NaN;
                    if n>=1200,at1200=h(1200);centerAt1200=g(1200);end
                    valid=logical(info.converged && info.AC_reconRelative<=cfg.maxAcLoss);
                    ctol=max(cfg.fs*cfg.tol,.001);
                    passedBy1200=logical(info.converged && n<=1200);
                    rr=table(string(name),string(labels{wi}),string(cfg.methods{k}),K,alpha,...
                        size(X,1),n,logical(info.converged),passedBy1200,valid,...
                        info.delta,info.centerDeltaHz,at1200,centerAt1200,...
                        info.AC_reconRelative,info.spectralResidual,secs,ctol,...
                        'VariableNames',{'fileName','windowRule','method','K','alpha','samples',...
                        'iterations','converged','convergedBy1200','solverValidAt2400',...
                        'finalDelta','finalCenterDeltaHz','deltaAt1200','centerDeltaAt1200Hz',...
                        'ACreconRelative','spectralResidual','seconds','centerToleranceHz'});
                    results=appendRows(results,rr);
                    hh=unique([cfg.logIterations(cfg.logIterations<=n) n]);
                    for j=1:numel(hh)
                        ix=hh(j);
                        histRow=table(string(name),string(labels{wi}),string(cfg.methods{k}),...
                            K,alpha,ix,h(ix),g(ix),info.residualHistory(ix),...
                            logical(h(ix)<cfg.tol && g(ix)<ctol),...
                            'VariableNames',{'fileName','windowRule','method','K','alpha',...
                            'iteration','delta','centerDeltaHz','spectralResidual','instantaneousTolerancePass'});
                        history=appendRows(history,histRow);
                    end
                    fprintf('    -> iter=%d converged=%d delta=%.3g center=%.4gHz AC=%.3f%%\n',...
                        n,info.converged,info.delta,info.centerDeltaHz,100*info.AC_reconRelative);
                    clear modes cent;
                catch ME
                    issues=appendRows(issues,newIssue(name,labels{wi},ME.message));
                    warning('E08F:Solver','%s %s: %s',name,cfg.methods{k},ME.message);
                end
                saveTables();
            end
        end
    catch ME
        issues=appendRows(issues,newIssue(name,'file_read',ME.message));
        warning('E08F:File','%s: %s',name,ME.message);
    end
    saveTables();
end
save(fullfile(cfg.outputDir,'E08F_config_frozen.mat'),'cfg');
fid=fopen(fullfile(cfg.outputDir,'E08F_说明.txt'),'w','n','UTF-8');
assert(fid>=0,'Cannot open summary');
fprintf(fid,'E08F 求解器诊断，不是FVIM寻优结果，也不是新的独立盲测。\n');
fprintf(fid,'固定3组E08E参数，保留tau=0、tol=1e-5、AC阈值0.15；2400迭代仅用于诊断收敛轨迹。\n');
fprintf(fid,'W9-60和Z9-3已用于算法问题诊断，不得再宣称全新盲测。\n');
fprintf(fid,'旧窗口为60%%至95%%且最多3秒；新诊断窗口为前0.5秒后至95%%时刻范围的最后最多3秒。\n');
fprintf(fid,'若更改正式求解迭代预算/窗口，所有方法必须同步更改、重新锁定训练配置并重新寻优，再用另一组未看过的数据独立验证。\n');
fprintf(fid,'仅凭本次求解器是否收敛，不判断角速度RMSE。\n');
fclose(fid);
fprintf('\nE08F finished: %s\n',cfg.outputDir);

    function saveTables()
        if ~isempty(windows),writetable(windows,fullfile(cfg.outputDir,'E08F_windows.csv'));end
        if ~isempty(results),writetable(results,fullfile(cfg.outputDir,'E08F_solver.csv'));end
        if ~isempty(history),writetable(history,fullfile(cfg.outputDir,'E08F_trajectory.csv'));end
        if ~isempty(issues),writetable(issues,fullfile(cfg.outputDir,'E08F_issues.csv'));end
    end
end

function p=locateOne(root,name)
a=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'));
    dir(fullfile(root,'*.TXT'));dir(fullfile(root,'*.txt'))];
match=a(arrayfun(@(x)strcmpi(x.name,name),a));
if ~isempty(match)
    paths=string(fullfile({match.folder},{match.name}));
    [~,u]=unique(lower(paths),'stable');match=match(u);
end
assert(numel(match)==1,'Need exactly one %s, found %d',name,numel(match));
p=fullfile(match.folder,match.name);
end

function T=appendRows(T,row)
if isempty(T),T=row;else,T=[T;row];end
end
function row=newIssue(name,rule,msg)
row=table(string(name),string(rule),string(msg),...
    'VariableNames',{'fileName','windowRule','message'});
end
