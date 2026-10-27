function E08H_run()
% E08-H predeclared K-reduction diagnostic, NOT FVIM and NOT blind test.
% Stage 1: 2 training files x 12 (K,alpha) combos; stage 2: 2 known dev files x selected + fixed references.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
cfg=E08H_config(); assert(isfolder(cfg.dataRoot),'请修改E08H_config.m中的cfg.dataRoot');
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
result=table();issues=table();selected=table();
for fi=1:numel(cfg.trainFiles)
 filename=cfg.trainFiles{fi};
 try
  [X,t,freqs]=readWindow(cfg,filename);
 catch ME
  issues=append(issues,issue(filename,'train_read',ME.message));saveProgress();continue
 end
 fprintf('\nTRAIN %s (%d samples)\n',filename,size(X,1));
 for ai=1:numel(cfg.alphaValues)
  for ki=1:numel(cfg.Kvalues)
   K=cfg.Kvalues(ki);alpha=cfg.alphaValues(ai);
   row=evaluate(cfg,X,t,freqs,filename,'train',K,alpha);
   result=append(result,row);saveProgress();
   fprintf('  K%d alpha=%.3f strict=%d ampPhase=%d score=%g\n',K,alpha,row.solverValid,row.fidelityPass,row.score);
  end
 end
end
if ~isempty(result)
 train=result(strcmp(result.stage,'train'),:);
 if height(train)>0
  for ki=1:numel(cfg.Kvalues)
   for ai=1:numel(cfg.alphaValues)
    K=cfg.Kvalues(ki);alpha=cfg.alphaValues(ai);
    p=train(train.K==K & abs(train.alpha-alpha)<1e-7,:);
    if height(p)~=numel(cfg.trainFiles),continue;end
    feasible=all(p.solverValid & p.fidelityPass) && all(isfinite(p.score));
    s=Inf;if feasible,s=mean(p.score);end
    sel=table(K,alpha,logical(feasible),s,...
     min(p.minRetentionPct),max(p.maxPhaseDeg),max(p.maxRelativePhaseDeg),...
     'VariableNames',{'K','alpha','feasibleBothTrainingFiles','meanScore','minRetentionPct','maxPhaseDeg','maxRelativePhaseDeg'});
    selected=append(selected,sel);
   end
  end
 end
end
best=zeros(0,2);
if ~isempty(selected)
 valid=selected(selected.feasibleBothTrainingFiles,:);
 if ~isempty(valid)
  valid=sortrows(valid,{'meanScore','K','alpha'});
  best=[valid.K(1) valid.alpha(1)];
 end
end
if isempty(best)
 fprintf('\nNO eligible candidate: training demands strict solver convergence and 99..101%% / 0.2deg conservation.\n');
 fprintf('This is a negative result, NOT a reason to tune on W9-60 or Z9-3.\n');
else
 fprintf('\nLOCKED from TRAIN: K=%d alpha=%.6f\n',best(1),best(2));
end
% Check development files after training-only selection. These were already seen in E08G.
C=unique([cfg.referenceCandidates;best],'rows','stable');
for fi=1:numel(cfg.developmentFiles)
 filename=cfg.developmentFiles{fi};
 try
  [X,t,freqs]=readWindow(cfg,filename);
 catch ME
  issues=append(issues,issue(filename,'development_read',ME.message));saveProgress();continue
 end
 fprintf('\nDEVELOPMENT %s (%d samples, NOT blind)\n',filename,size(X,1));
 for j=1:size(C,1)
  row=evaluate(cfg,X,t,freqs,filename,'development',C(j,1),C(j,2));
  result=append(result,row);saveProgress();
  fprintf('  K%d alpha=%.3f strict=%d fidelity=%d score=%g\n',C(j,1),C(j,2),row.solverValid,row.fidelityPass,row.score);
 end
end
save(fullfile(cfg.outputDir,'E08H_config_frozen.mat'),'cfg','best');
if isempty(best),verdict='NO_TRAINING_FEASIBLE';else,verdict='TRAINING_CANDIDATE_SELECTED';end
fid=fopen(fullfile(cfg.outputDir,'E08H_decision.txt'),'w','n','UTF-8');
if fid>=0
 fprintf(fid,'E08-H: 模态数降阶严格收敛验证。非自动FVIM，不是盲测。\n');
 fprintf(fid,'原始求解器eta=1, tol=1e-5, maxIter=2400，AC<=15%%。训练文件W9-10与Z9-1。\n');
 fprintf(fid,'规则：训练两文件都必须严格收敛、重构合格、载波/边带幅值99..101%%及相位<=0.2度。\n');
 fprintf(fid,'选参结果：%s\n',verdict);
 if ~isempty(best),fprintf(fid,'训练选出 K=%d alpha=%.9f，不依赖开发集W9-60/Z9-3。\n',best(1),best(2));end
 fprintf(fid,'W9-60和Z9-3均已用于前序研究，只能做开发复核，不得称全新独立盲测。\n');
 fprintf(fid,'未收敛候选的score设为Inf，不能纳入正式角速度RMSE。\n');
 fclose(fid);
end
saveProgress();fprintf('\nE08H completed. %s\n',cfg.outputDir);
 function saveProgress()
  if ~isempty(result),writetable(result,fullfile(cfg.outputDir,'E08H_cases.csv'));end
  if ~isempty(selected),writetable(selected,fullfile(cfg.outputDir,'E08H_training_selection.csv'));end
  if ~isempty(issues),writetable(issues,fullfile(cfg.outputDir,'E08H_issues.csv'));end
 end
end
function row=evaluate(cfg,X,t,freqs,filename,stage,K,alpha)
ticRun=tic;
conv=false;valid=false;pass=false;it=NaN;delta=NaN;cgap=NaN;ac=NaN;
sc=Inf;minRet=NaN;maxPhase=NaN;maxRel=NaN;noiseA=NaN;noiseS=NaN;
domMode=NaN;domEnergy=NaN;domChange=NaN;carrierEnergy=NaN;
errMsg="";
try
 [U,cent,info,snaps]=E08G_mvmd_damped(X,cfg.fs,K,alpha,cfg.maxIter,cfg.tol,cfg.tau,cfg.rollHz,1,cfg.maxIter);
 conv=logical(info.converged);it=info.iterations;delta=info.delta;
 cgap=info.centerDeltaHz;ac=info.ACreconRelative;
 valid=conv && ac<=cfg.maxAcLoss;
 sn=snaps(end);[domChange,j]=max(sn.modeGlobalDeltaShare);
 domMode=j;domEnergy=sn.modeEnergyShare(j);
 [~,j9]=min(abs(sn.centerHz-cfg.rollHz));carrierEnergy=sn.modeEnergyShare(j9);
 if valid
  [Y,~]=gyro_reconstruct(U,cent,X,cfg,'paper');
  [score,diag]=E07D_score(X,Y,U,t,freqs,freqs,cfg);
  mm=diag.metrics;minRet=mm.minSidebandRetentionPct;
  maxPhase=mm.maxPhaseShiftDeg;maxRel=mm.maxRelativePhaseShiftDeg;
  noiseA=mm.nonTargetSuppressionDb(1);noiseS=mm.nonTargetSuppressionDb(2);
  pass=logical(mm.sidebandPass)&&mm.maxRelativePhaseShiftDeg<=cfg.pass.maxSidebandPhaseShiftDeg;
  if pass,sc=score;end
 end
catch ME
 errMsg=string(ME.message);
end
row=table(string(filename),string(stage),K,alpha,logical(conv),logical(valid),logical(pass),...
 it,delta,cgap,ac,sc,minRet,maxPhase,maxRel,noiseA,noiseS,domMode,domEnergy,domChange,carrierEnergy,...
 toc(ticRun),errMsg,'VariableNames',{'fileName','stage','K','alpha','converged','solverValid',...
 'fidelityPass','iterations','finalModeDelta','centerGapHz','ACreconRelative','score','minRetentionPct',...
 'maxPhaseDeg','maxRelativePhaseDeg','proxyA_dB','proxyS_dB','dominantChangeModeIndex',...
 'dominantChangeModeEnergyFraction','dominantChangeFraction','carrierModeEnergyFraction','seconds','errorMessage'});
end
function [X,t,freqs]=readWindow(cfg,filename)
p=locateOne(cfg.dataRoot,filename);r=read_gyro_txt(p);
assert(r.invalidRows==0 && r.headerBlockRuleMatches,'TXT data format failed');
N=size(r.data,1);stop=min(N,floor(cfg.endFraction*N));
start=max(1+round(cfg.skipStartSec*cfg.fs),stop-cfg.lastSeconds*cfg.fs+1);
X=double(r.data(start:stop,[4 7]));
assert(size(X,1)>=cfg.minSamples,sprintf('%s too short after window selection',filename));
t=(0:size(X,1)-1)'/cfg.fs;
tok=regexp(upper(filename),'^Z(\d+)-(\d+)\.TXT$','tokens','once');
if ~isempty(tok)
  modulationHz=str2double(tok{2});
  freqs=unique([cfg.rollHz-modulationHz,cfg.rollHz,cfg.rollHz+modulationHz]);
  assert(all(freqs>0 & freqs<cfg.fs/2),'Invalid modulation sidebands');
else
  freqs=cfg.rollHz;
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
function t=issue(file,stage,msg)
t=table(string(file),string(stage),string(msg),...
 'VariableNames',{'fileName','stage','message'});
end
function T=append(T,r)
if isempty(T),T=r;else,T=[T;r];end
end
