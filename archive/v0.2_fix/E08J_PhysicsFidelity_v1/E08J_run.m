function E08J_run()
% E08-J: small, controlled solver and fidelity-rule experiment.
% NOT a new FVIM or a blind validation. Only training data choose candidates.
cfg=E08J_config();assert(isfolder(cfg.dataRoot),'修改 E08J_config.m 的数据目录');
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
save(fullfile(cfg.outputDir,'E08J_frozen_config.mat'),'cfg');
allCases=table();allBands=table();allModes=table();errors=table();
files=cfg.trainFiles;
if cfg.enableKnownDevelopment,files=[files cfg.developmentFiles];end
for fi=1:numel(files)
 filename=files{fi};stage='training';
 if fi>numel(cfg.trainFiles),stage='known_development';end
 try
  [X,t,freqs]=getWindow(cfg,filename);
 catch ME
  errors=appendRow(errors,table(string(filename),string(stage),NaN,NaN,string(ME.message),...
   'VariableNames',{'fileName','stage','K','alpha','errorMessage'}));flush();continue;
 end
 fprintf('E08J %s (%s): %d samples; model frequencies: %s Hz\n',...
  filename,stage,numel(t),mat2str(freqs));
 for ia=1:numel(cfg.alphas)
  for ik=1:numel(cfg.Ks)
   K=cfg.Ks(ik);alpha=cfg.alphas(ia);
   try
    ticStart=tic;
    [modes,cent,info,~]=E08G_mvmd_damped(X,cfg.fs,K,alpha,...
      cfg.maxIter,cfg.tol,cfg.tau,cfg.rollHz,cfg.eta,cfg.maxIter);
    solverTime=toc(ticStart);
    solverValid=logical(info.converged)&&(info.ACreconRelative<=cfg.maxAcLoss);
    if solverValid
      [hard,imf]=gyro_reconstruct(modes,cent,X,cfg,'paper');
      full=sum(modes,3)+mean(X,1);
      representations={'hard_paper',hard;'all_modes',full};
      for v=1:size(representations,1)
       reconstruction=representations{v,1};Y=representations{v,2};
       [qual,b]=E08J_physics_metrics(X,Y,t,freqs,cfg);
       first=table(string(filename),string(stage),K,alpha,string(reconstruction),...
          logical(solverValid),logical(info.converged),info.iterations,info.delta,...
          info.centerDeltaHz,info.ACreconRelative,solverTime,logical(qual.oldPass),...
          logical(qual.newPass),qual.minimumStrongRetention,...
          qual.maximumStrongRetention,qual.maximumWeakComplexFraction,...
          qual.maximumStrongPhaseDeg,qual.maximumStrongRelativePhaseDeg,...
          qual.minimumProxyDb,qual.meanProxyDb,...
          'VariableNames',{'fileName','stage','K','alpha','reconstruction','solverValid',...
           'converged','iterations','modeDelta','centerGapHz','ACreconRelative',...
           'seconds','oldPhysicsPass','newPhysicsPass','minStrongRetentionPct',...
           'maxStrongRetentionPct','maxWeakComplexFraction','maxStrongPhaseDeg',...
           'maxStrongRelativePhaseDeg','minProxyDb','meanProxyDb'});
       allCases=appendRow(allCases,first);
       b.fileName=repmat(string(filename),height(b),1);
       b.stage=repmat(string(stage),height(b),1);
       b.K=repmat(K,height(b),1);b.alpha=repmat(alpha,height(b),1);
       b.reconstruction=repmat(string(reconstruction),height(b),1);
       allBands=appendRow(allBands,b);
      end
      m=table(repmat(string(filename),height(imf),1),repmat(K,height(imf),1),...
       repmat(alpha,height(imf),1),imf.mode,imf.center_Hz,imf.energy_share,...
       imf.cosine_abs,imf.sampEn,imf.class_0drop_1mix_2keep,...
       'VariableNames',{'fileName','K','alpha','IMF','centerHz','energyShare','CS','SampEn','classification'});
      allModes=appendRow(allModes,m);
    else
      row=table(string(filename),string(stage),K,alpha,"not_reconstructed",...
       logical(solverValid),logical(info.converged),info.iterations,info.delta,...
       info.centerDeltaHz,info.ACreconRelative,solverTime,false,false,...
       NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
       'VariableNames',{'fileName','stage','K','alpha','reconstruction','solverValid',...
        'converged','iterations','modeDelta','centerGapHz','ACreconRelative',...
        'seconds','oldPhysicsPass','newPhysicsPass','minStrongRetentionPct',...
        'maxStrongRetentionPct','maxWeakComplexFraction','maxStrongPhaseDeg',...
        'maxStrongRelativePhaseDeg','minProxyDb','meanProxyDb'});
      allCases=appendRow(allCases,row);
    end
    fprintf(' K=%d alpha=%.4f converged=%d iter=%d old/new=%s\n',...
       K,alpha,info.converged,info.iterations,showFidelity(allCases,filename,K,alpha));
    flush();
   catch ME
    fprintf(' ERROR %s K=%d a=%.3f: %s\n',filename,K,alpha,ME.message);
    errors=appendRow(errors,table(string(filename),string(stage),K,alpha,string(ME.message),...
     'VariableNames',{'fileName','stage','K','alpha','errorMessage'}));flush();
   end
  end
 end
end
% Never use known_development results to select parameters.
selection=table();
if ~isempty(allCases)
 M=allCases(allCases.stage=="training" & allCases.reconstruction=="hard_paper",:);
 for ia=1:numel(cfg.alphas)
  for ik=1:numel(cfg.Ks)
   K=cfg.Ks(ik);alpha=cfg.alphas(ia);
   sub=M(M.K==K & abs(M.alpha-alpha)<1e-6,:);
   allSeen=height(sub)==numel(cfg.trainFiles);
   feasible=allSeen && all(sub.solverValid) && all(sub.newPhysicsPass);
   robustProxy=NaN;meanProxy=NaN;score=-Inf;
   if feasible
     robustProxy=min(sub.minProxyDb);meanProxy=mean(sub.meanProxyDb);
     score=robustProxy+0.1*meanProxy-0.001*K;
   end
   selection=appendRow(selection,table(K,alpha,allSeen,feasible,robustProxy,...
     meanProxy,score,'VariableNames',{'K','alpha','allTrainingSeen',...
     'newPhysicsFeasible','worstChannelProxyDb','meanProxyDb','exploratoryScore'}));
  end
 end
end
if ~isempty(selection),writetable(selection,fullfile(cfg.outputDir,'E08J_train_selection.csv'));end
fid=fopen(fullfile(cfg.outputDir,'E08J_decision.txt'),'w','n','UTF-8');
if fid>=0
 fprintf(fid,'E08-J 仅改物理评价规则；原始 MVMD / 重构完全不变。\n');
 fprintf(fid,'旧规则：所有频点幅值99~101%%、显著相位<=0.2deg及双路相对相位。\n');
 fprintf(fid,'新规则：强分量(>=同通道最强目标分量10%%)仍需幅值1%%和相位0.2deg；弱分量用复系数绝对误差<=最强幅值1%%；双通道强/强才检验相对相位0.2deg。\n');
 fprintf(fid,'此10%%划分在已见Z9-1基础上制定，是开发阶段探索，不能称盲测确认。\n');
 fprintf(fid,'求解器tol1e-5、2400次、AC<=15%%不变，未收敛不可正式比较。\n');
 fprintf(fid,'训练仅W9-10、Z9-1；W9-60/Z9-3仅已知开发复核，禁止反向选参。\n');
 fprintf(fid,'评分为非目标残差代理，非真实SNR；这里不输出编码器角速度RMSE/姿态RMSE。\n');
 if isempty(selection)||~any(selection.newPhysicsFeasible)
  fprintf(fid,'训练结论：未选出通过全部保真和收敛条件的候选；不得假定FVIM自动寻优成功。\n');
 else
  feasible=selection(selection.newPhysicsFeasible,:);
  [~,best]=max(feasible.exploratoryScore);pick=feasible(best,:);
  fprintf(fid,'开发阶段候选 K=%d alpha=%.9f;仅为固定候选比较，不是FVIM自动寻优。\n',pick.K,pick.alpha);
 end
 fclose(fid);
end
flush();fprintf('E08J finished: %s\n',cfg.outputDir);
 function flush()
  if ~isempty(allCases),writetable(allCases,fullfile(cfg.outputDir,'E08J_cases.csv'));end
  if ~isempty(allBands),writetable(allBands,fullfile(cfg.outputDir,'E08J_band_detail.csv'));end
  if ~isempty(allModes),writetable(allModes,fullfile(cfg.outputDir,'E08J_modes.csv'));end
  if ~isempty(errors),writetable(errors,fullfile(cfg.outputDir,'E08J_issues.csv'));end
 end
end
function str=showFidelity(tbl,filename,K,alpha)
 match=tbl(tbl.fileName==string(filename)&tbl.K==K&abs(tbl.alpha-alpha)<1e-6&tbl.reconstruction=="hard_paper",:);
 if isempty(match),str='no valid reconstruction';else
  str=sprintf('%d / %d',match.oldPhysicsPass(end),match.newPhysicsPass(end));end
end
function [X,t,freqs]=getWindow(cfg,filename)
 files=[dir(fullfile(cfg.dataRoot,'**','*.TXT'));dir(fullfile(cfg.dataRoot,'**','*.txt'));...
 dir(fullfile(cfg.dataRoot,'*.TXT'));dir(fullfile(cfg.dataRoot,'*.txt'))];
 hit=files(arrayfun(@(d)strcmpi(d.name,filename),files));
 if ~isempty(hit)
  paths=string(fullfile({hit.folder},{hit.name}));[~,ia]=unique(lower(paths),'stable');hit=hit(ia);
 end
 assert(numel(hit)==1,'数据文件 %s 找到%d个，必须唯一',filename,numel(hit));
 dat=read_gyro_txt(fullfile(hit.folder,hit.name));
 assert(dat.invalidRows==0 && dat.headerBlockRuleMatches,'TXT数据校验失败');
 N=size(dat.data,1);last=min(N,floor(cfg.endFraction*N));
 first=max(1+round(cfg.skipStartSec*cfg.fs),last-cfg.lastSeconds*cfg.fs+1);
 assert(last-first+1>=cfg.minSamples,'有效窗口不足要求');
 X=double(dat.data(first:last,[4 7]));t=(0:size(X,1)-1)'/cfg.fs;
 tok=regexp(upper(filename),'^Z(\d+)-(\d+)\.TXT$','tokens','once');
 if isempty(tok),freqs=cfg.rollHz;else
  dHz=str2double(tok{2});freqs=unique([cfg.rollHz-dHz,cfg.rollHz,cfg.rollHz+dHz]);
 end
 assert(all(freqs>0 & freqs<cfg.fs/2),'边带频率越界');
end
function T=appendRow(T,row)
 if isempty(T),T=row;else,T=[T;row];end
end
