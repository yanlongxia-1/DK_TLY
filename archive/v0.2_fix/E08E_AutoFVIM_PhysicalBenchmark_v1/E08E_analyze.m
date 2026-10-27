function E08E_analyze()
%E08E_ANALYZE 只读真实输出，在同一W文件上的有效算法之间计算RMSE变化；不参与寻优。
cfg=E08E_config();p=fullfile(cfg.outputDir,'E08E_test_metrics.csv');
assert(isfile(p),'先运行E08E_run，缺少:%s',p);
T=readtable(p,'TextType','string');
fileList=unique(T.fileName,'stable');summary=table();
for i=1:numel(fileList)
    fn=fileList(i);a=T(T.fileName==fn,:);
    get=@(n) a(a.method==string(n),:);
    fixed=get('fixed_baseline');ent=get('entropy_auto');new=get('physics_auto');
    assert(height(fixed)==1&&height(ent)==1&&height(new)==1,'同文件缺少一种方法行');
    if startsWith(upper(fn),'W')
        isComparable=logical(fixed.solverValid && new.solverValid && new.fidelityPass && ...
            isfinite(fixed.RMSE_deg_s) && isfinite(new.RMSE_deg_s));
        deltaFixed=NaN;if isComparable
            deltaFixed=100*(fixed.RMSE_deg_s-new.RMSE_deg_s)/max(fixed.RMSE_deg_s,eps);
        end
        entComparable=logical(ent.solverValid && new.solverValid && new.fidelityPass && ...
            isfinite(ent.RMSE_deg_s) && isfinite(new.RMSE_deg_s));
        deltaEntropy=NaN;if entComparable
            deltaEntropy=100*(ent.RMSE_deg_s-new.RMSE_deg_s)/max(ent.RMSE_deg_s,eps);
        end
        fprintf('%s: frozenValid=%d entropyValid=%d physicsValid=%d | RMSE改善vs冻结=%.4g%%, vs熵=%.4g%%\n',...
            fn,fixed.solverValid,ent.solverValid,new.solverValid,deltaFixed,deltaEntropy);
    else
        isComparable=false;entComparable=false;deltaFixed=NaN;deltaEntropy=NaN;
        fprintf('%s: Z动态工况，仅评价求解有效率/边带保真，暂不计算恒速角速度RMSE\n',fn);
    end
    row=table(fn,logical(fixed.solverValid),logical(ent.solverValid),logical(new.solverValid),...
        logical(new.fidelityPass),isComparable,deltaFixed,entComparable,deltaEntropy,...
        'VariableNames',{'fileName','fixedValid','entropyValid','physicsValid','physicsFidelityPass',...
        'pairedWithFixed','RMSE_improvementVsFixed_pct','pairedWithEntropy','RMSE_improvementVsEntropy_pct'});
    if isempty(summary),summary=row;else,summary=[summary;row];end %#ok<AGROW>
end
writetable(summary,fullfile(cfg.outputDir,'E08E_paired_summary.csv'));
fprintf('只读取现有测试指标，不会重新寻优，输出E08E_paired_summary.csv。\n');
end
