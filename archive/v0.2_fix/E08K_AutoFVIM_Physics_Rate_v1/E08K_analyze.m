function E08K_analyze()
% Post-run analysis only: no optimization, no test-driven parameter tuning.
cfg=E08K_config();path=fullfile(cfg.outputDir,'E08K_test_metrics.csv');
assert(isfile(path),'请先运行 E08K_run');
T=readtable(path,'TextType','string');
fid=fopen(fullfile(cfg.outputDir,'E08K_comparison.txt'),'w','n','UTF-8');
assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'E08-K: W9-60与Z9-3历史已用于开发，不是新盲测。\n');
fprintf(fid,'新物理保真通过不等于真值SNR或姿态RMSE提升。\n');
for ii=1:height(T)
 fprintf(fid,'%s %s valid=%d fidelity=%d rmseDegS=%.8g\n',...
  T.fileName(ii),T.method(ii),T.solverValid(ii),T.fidelityPass(ii),T.RMSE_deg_s(ii));
end
w=T(contains(upper(T.fileName),'W9-60') & T.solverValid & T.fidelityPass & isfinite(T.RMSE_deg_s),:);
if ~isempty(w)
 raw=w(w.method=="raw",:);
 if ~isempty(raw)
  for i=1:height(w)
   fprintf(fid,'Compared to RAW, %s deltaRMSEpct=%.4f\n',w.method(i),...
     100*(raw.RMSE_deg_s(1)-w.RMSE_deg_s(i))/raw.RMSE_deg_s(1));
  end
 end
end
fprintf('E08-K summary saved: %s\n',fullfile(cfg.outputDir,'E08K_comparison.txt'));
end
