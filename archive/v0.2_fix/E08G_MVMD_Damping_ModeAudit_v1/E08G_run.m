function E08G_run()
% E08-G: isolate slow/nonconvergent mode dynamics. W9-60/Z9-3 DEVELOPMENT ONLY.
% No K/alpha search, no encoder RMSE, no tuning from previously held-out W9-60.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'lib'));
cfg=E08G_config();
assert(isfolder(cfg.dataRoot),'Set cfg.dataRoot in E08G_config.m');
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
solverRows=table();modeRows=table();checkpointRows=table();stabilityRows=table();issueRows=table();
for fi=1:numel(cfg.filenames)
 filename=cfg.filenames{fi};
 try
  path=locateOne(cfg.dataRoot,filename);
  r=read_gyro_txt(path);
  assert(r.invalidRows==0 && r.headerBlockRuleMatches,'Malformed TXT file');
  N=size(r.data,1);stop=min(N,floor(.95*N));
  start=max(1+round(cfg.skipStartSec*cfg.fs),stop-cfg.lastSeconds*cfg.fs+1);
  X=double(r.data(start:stop,[4 7]));
  assert(size(X,1)>=2000,'Insufficient valid samples after tail-window selection');
  t=(0:size(X,1)-1)'/cfg.fs;
  fprintf('\nE08G %s rows=%d window %d:%d\n',filename,N,start,stop);
  for k=1:size(cfg.candidates,1)
   K=cfg.candidates(k,1);alpha=cfg.candidates(k,2);method=cfg.names{k};
   for e=1:numel(cfg.relaxValues)
    eta=cfg.relaxValues(e);tRun=tic;
    try
     [~,~,info,snaps]=E08G_mvmd_damped(X,cfg.fs,K,alpha,...
       cfg.maxIter,cfg.tol,cfg.tau,cfg.rollHz,eta,cfg.checkpoints);
     secs=toc(tRun);
     valid=info.converged && info.ACreconRelative<=cfg.maxAcLoss;
     sr=table(string(filename),string(method),K,alpha,eta,...
       size(X,1),info.iterations,logical(info.converged),logical(valid),...
       info.delta,info.centerDeltaHz,info.centerStepHz,info.ACreconRelative,secs,...
       'VariableNames',{'fileName','method','K','alpha','centerRelaxation',...
       'samples','iterations','converged','solverValid',...
       'finalModeDelta','finalCenterGapHz','finalCenterStepHz','ACreconRelative','seconds'});
     solverRows=addRows(solverRows,sr);
     yy=cell(numel(snaps),1); hh=zeros(numel(snaps),2);cc=cell(numel(snaps),1);
     for s=1:numel(snaps)
      a=snaps(s); [y,classes]=gyro_reconstruct(a.modes,a.centerHz,X,cfg,'paper');
      yy{s}=y;cc{s}=classes.class_0drop_1mix_2keep;
      h=zeros(1,2);
      for ch=1:2
       x=y(:,ch)-mean(y(:,ch));
       h(ch)=2/numel(t)*sum(x.*exp(-1i*2*pi*cfg.rollHz*t));
      end
      hh(s,:)=h;
      kept=sum(cc{s}==2);mixed=sum(cc{s}==1);
      cp=table(string(filename),string(method),K,alpha,eta,a.iteration,a.delta,...
        a.centerDeltaHz,a.centerStepHz,a.ACreconRelative,kept,mixed,...
        abs(h(1)),abs(h(2)),wrapDeg(angle(h(2)/h(1))*180/pi),...
        'VariableNames',{'fileName','method','K','alpha','centerRelaxation',...
        'iteration','modeDelta','centerGapHz','centerStepHz','ACreconRelative',...
        'keptModes','mixedModes','carrierAmpA','carrierAmpS','relativePhaseDeg'});
      checkpointRows=addRows(checkpointRows,cp);
      for j=1:K
       mr=table(string(filename),string(method),K,alpha,eta,a.iteration,j,...
         a.centerHz(j),a.modeEnergyShare(j),a.modeLocalDelta(j),a.modeGlobalDeltaShare(j),cc{s}(j),...
         'VariableNames',{'fileName','method','K','alpha','centerRelaxation','iteration','modeIndex',...
         'centerHz','energyFraction','localModeDelta','globalChangeFraction','reconstructionClass'});
       modeRows=addRows(modeRows,mr);
      end
     end
     % Final-vs-earlier output changes are DIAGNOSTICS, not proof of convergence.
     last=numel(snaps); yLast=yy{last};hLast=hh(last,:);
     for s=1:last-1
      for ch=1:2
       sig=yy{s}(:,ch);baseline=yLast(:,ch);
       changeRms=norm(sig-baseline)/max(norm(baseline-mean(baseline)),eps);
       ampDiff=100*(abs(hh(s,ch))/max(abs(hLast(ch)),eps)-1);
       phaseDiff=wrapDeg((angle(hh(s,ch))-angle(hLast(ch)))*180/pi);
       classFlips=sum(cc{s}~=cc{last});
       st=table(string(filename),string(method),eta,snaps(s).iteration,...
         snaps(last).iteration,ch,changeRms,ampDiff,phaseDiff,classFlips,...
         'VariableNames',{'fileName','method','centerRelaxation','earlyIter',...
         'referenceIter','channel','relativeRmsDifference','carrierAmplitudeChangePct',...
         'carrierPhaseChangeDeg','classChanges'});
       stabilityRows=addRows(stabilityRows,st);
      end
     end
     fprintf('  %-14s eta=%.2f iter=%d pass=%d delta=%.3g center=%.4gHz AC=%.3f%%\n',...
        method,eta,info.iterations,valid,info.delta,info.centerDeltaHz,100*info.ACreconRelative);
    catch ME
     issueRows=addRows(issueRows,table(string(filename),string(method),eta,string(ME.message),...
        'VariableNames',{'fileName','method','centerRelaxation','message'}));
     warning('E08G:Failure','%s / %s / %.2f: %s',filename,method,eta,ME.message);
    end
    persist();
   end
  end
 catch ME
  issueRows=addRows(issueRows,table(string(filename),string('file_read'),NaN,string(ME.message),...
      'VariableNames',{'fileName','method','centerRelaxation','message'}));
  warning('E08G:File','%s: %s',filename,ME.message);
 end
 persist();
end
save(fullfile(cfg.outputDir,'E08G_config_frozen.mat'),'cfg');
fid=fopen(fullfile(cfg.outputDir,'E08G_decision_notes.txt'),'w','n','UTF-8');
if fid>=0
 fprintf(fid,'E08G是开发阶段求解器消融：原中心更新eta=1与预定阻尼eta=0.5；不是FVIM优化或盲测。\n');
 fprintf(fid,'严格tol=1e-5、中心0.01Hz、连续3次、AC阈值0.15保持不变。\n');
 fprintf(fid,'非收敛迭代快照仅用于机制诊断，不能用于正式角速度RMSE比较。\n');
 fprintf(fid,'若阻尼改善，必须在训练/开发集确定eta并冻结，对未使用文件进行统一预算检验。\n');
 fclose(fid);
end
fprintf('E08G completed: %s\n',cfg.outputDir);

 function persist()
  if ~isempty(solverRows),writetable(solverRows,fullfile(cfg.outputDir,'E08G_solver.csv'));end
  if ~isempty(modeRows),writetable(modeRows,fullfile(cfg.outputDir,'E08G_modes.csv'));end
  if ~isempty(checkpointRows),writetable(checkpointRows,fullfile(cfg.outputDir,'E08G_checkpoints.csv'));end
  if ~isempty(stabilityRows),writetable(stabilityRows,fullfile(cfg.outputDir,'E08G_output_stability.csv'));end
  if ~isempty(issueRows),writetable(issueRows,fullfile(cfg.outputDir,'E08G_issues.csv'));end
 end
end
function v=wrapDeg(v),v=mod(v+180,360)-180;end
function T=addRows(T,r)
if isempty(T)
 T=r;
else
 T=[T;r];
end
end
function p=locateOne(root,name)
a=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'));...
 dir(fullfile(root,'*.TXT'));dir(fullfile(root,'*.txt'))];
match=a(arrayfun(@(x)strcmpi(x.name,name),a));
if ~isempty(match)
 paths=string(fullfile({match.folder},{match.name}));
 [~,u]=unique(lower(paths),'stable');match=match(u);
end
assert(numel(match)==1,'Expected exactly one %s, found %d',name,numel(match));
p=fullfile(match.folder,match.name);
end
