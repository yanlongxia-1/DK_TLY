function E08I_run()
% 预登记目标：一、标记Z9-1超出101%的真实边带/通道；二、判定是MVMD分解还是硬重构造成。
% 三、检验低K在W9-60能否严格收敛；不是FVIM自动选参，不使用编码器真值。
cfg=E08I_config();
assert(isfolder(cfg.dataRoot), '请修改 E08I_config.m 中的 cfg.dataRoot');
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
stats=table(); bands=table(); modes=table(); issues=table();
items={ 'Z9-1.TXT','train_diagnostic',cfg.Zcandidates; ...
        'W9-60.TXT','known_development',cfg.Wcandidates};
for ii=1:size(items,1)
 file=items{ii,1}; stage=items{ii,2}; cand=items{ii,3};
 try
  [X,t,freqs]=readWindow(cfg,file);
 catch ME
  issues=addrow(issues,table(string(file),string(stage),string(ME.message), ...
   'VariableNames',{'fileName','stage','errorMessage'})); saveCsv(); continue;
 end
 fprintf('\nE08I %s %d点 / %s\n',file,size(X,1),stage);
 for k=1:size(cand,1)
  K=cand(k,1); alpha=cand(k,2); t0=tic;
  try
   [U,c,inf,~]=E08G_mvmd_damped(X,cfg.fs,K,alpha,cfg.maxIter,cfg.tol,cfg.tau,cfg.rollHz,cfg.eta,cfg.maxIter);
   solverOK=logical(inf.converged) && inf.ACreconRelative<=cfg.maxAcLoss;
   fprintf('  K%d a%.3f iter=%d converged=%d AC=%.4f%%\n',K,alpha,inf.iterations,inf.converged,100*inf.ACreconRelative);
   stat=table(string(file),string(stage),K,alpha,logical(inf.converged),logical(solverOK), ...
    inf.iterations,inf.delta,inf.centerDeltaHz,inf.ACreconRelative,toc(t0), ...
    'VariableNames',{'fileName','stage','K','alpha','converged','solverValid', ...
    'iterations','finalModeDelta','centerGapHz','ACreconRelative','seconds'});
   stats=addrow(stats,stat); saveCsv();
   if ~solverOK,continue;end
   [paper,det]=gyro_reconstruct(U,c,X,cfg,'paper');
   classes=det.class_0drop_1mix_2keep;
   fullRecon=sum(U,3)+mean(X,1);
   mixedRaw=zeros(size(X)); for j=1:numel(classes)
    if classes(j)>0,mixedRaw=mixedRaw+U(:,:,j);end
   end
   mixedRaw=mixedRaw+mean(X,1);
   variations={'hard_paper',paper;'all_modes',fullRecon;'mixed_without_wavelet',mixedRaw};
   for j=1:size(det,1)
    row=table(string(file),K,alpha,j,det.center_Hz(j),det.energy_share(j), ...
     det.cosine_abs(j),det.sampEn(j),classes(j), ...
     'VariableNames',{'fileName','K','alpha','IMF','centerHz','energyShare','CS','SampEn','classification'});
    modes=addrow(modes,row);
   end
   for v=1:size(variations,1)
    name=variations{v,1};Y=variations{v,2};
    m=E07B_candidate_metrics(X,Y,t,freqs,freqs,cfg);
    % sidebandPass还需单独检查两路相对相位，而非只检查单路。
    pass=logical(m.sidebandPass) && m.maxRelativePhaseShiftDeg<=cfg.pass.maxSidebandPhaseShiftDeg;
    fprintf('    %s min=%.4f%% MAX=%.4f%% phase=%.4fdeg rel=%.4fdeg PASS=%d\n',...
     name,min(m.retentionPct(:)),max(m.retentionPct(:)),m.maxPhaseShiftDeg,m.maxRelativePhaseShiftDeg,pass);
    for ch=1:2
     for f=1:numel(freqs)
      % 使用与已有评价代码相同的联合频率最小二乘，统计原始幅值和拟合残差。
      [rawAmp,fitResid]=rawFitDiagnostics(X(:,ch),t,freqs,f,ch);
      badLow=m.retentionPct(ch,f)<cfg.pass.minRetentionPct;
      badHigh=m.retentionPct(ch,f)>cfg.pass.maxRetentionPct;
      badPhase=m.phaseValid(ch,f) && m.phaseAbsDeg(ch,f)>cfg.pass.maxSidebandPhaseShiftDeg;
      % 双通道相对相位变化是频点公共指标；不能把它当作独立通道误差。
      badRel=m.relativePhaseChangeDeg(f)>cfg.pass.maxSidebandPhaseShiftDeg;
      row=table(string(file),string(stage),K,alpha,string(name),ch,freqs(f),rawAmp,fitResid,...
       m.retentionPct(ch,f),m.phaseAbsDeg(ch,f),m.relativePhaseChangeDeg(f),...
       logical(m.phaseValid(ch,f)),logical(badLow),logical(badHigh),logical(badPhase),logical(badRel),...
       logical(pass),m.nonTargetSuppressionDb(ch),...
       'VariableNames',{'fileName','stage','K','alpha','reconstruction','channel','frequencyHz',...
       'rawFitAmplitude','rawFitResidualRMS','retentionPct','phaseChangeDeg','relativePhaseChangeDeg',...
       'phaseValid','lowAmpViolation','highAmpViolation','phaseViolation','relativePhaseViolation',...
       'allPhysicsPass','nonTargetSuppressionDb'});
      bands=addrow(bands,row);
     end
    end
   end
   saveCsv();
  catch ME
   issues=addrow(issues,table(string(file),string(stage),K,alpha,string(ME.message),...
    'VariableNames',{'fileName','stage','K','alpha','errorMessage'}));saveCsv();
   fprintf('  ERROR K%d alpha%.3f: %s\n',K,alpha,ME.message);
  end
 end
end
fid=fopen(fullfile(cfg.outputDir,'E08I_notes.txt'),'w','n','UTF-8');
if fid>=0
 fprintf(fid,'E08-I诊断：不能用W9-60选择K/alpha；所有阈值保持E08-H原值。\n');
 fprintf(fid,'Z9-1对照硬重构、全部模态重构、混叠模态不小波重构，找出幅值>101%%的频点。\n');
 fprintf(fid,'全模态通过、硬重构失败 => 优先排查分类/重构；全模态也失败 => 优先排查分解/重构残差。\n');
 fprintf(fid,'W9-60仅检查小K严格收敛，不将未收敛快照用于角速度RMSE。\n');
 fprintf(fid,'输出不包含真实编码器精度，也不是新的独立盲测。\n');
 fclose(fid);
end
saveCsv(); fprintf('\nE08I结束：%s\n',cfg.outputDir);
 function saveCsv()
  if ~isempty(stats),writetable(stats,fullfile(cfg.outputDir,'E08I_solver.csv'));end
  if ~isempty(bands),writetable(bands,fullfile(cfg.outputDir,'E08I_band_detail.csv'));end
  if ~isempty(modes),writetable(modes,fullfile(cfg.outputDir,'E08I_mode_classes.csv'));end
  if ~isempty(issues),writetable(issues,fullfile(cfg.outputDir,'E08I_issues.csv'));end
 end
end
function [rawAmp,residRMS]=rawFitDiagnostics(x,t,freqs,findex,~)
F=ones(numel(t),1);for q=1:numel(freqs)
 F=[F cos(2*pi*freqs(q)*t) sin(2*pi*freqs(q)*t)]; %#ok<AGROW>
end
b=F\x;rawAmp=hypot(b(2*findex),b(2*findex+1));
residRMS=sqrt(mean((x-F*b).^2));
end
function [X,t,freqs]=readWindow(cfg,filename)
p=locateOne(cfg.dataRoot,filename);r=read_gyro_txt(p);
assert(r.invalidRows==0 && r.headerBlockRuleMatches,'TXT读取格式校验失败');
N=size(r.data,1);last=min(N,floor(cfg.endFraction*N));
start=max(1+round(cfg.skipStartSec*cfg.fs),last-cfg.lastSeconds*cfg.fs+1);
X=double(r.data(start:last,[4 7]));
assert(size(X,1)>=cfg.minSamples,sprintf('%s 数据窗口不足',filename));
t=(0:size(X,1)-1)'/cfg.fs;
match=regexp(upper(filename),'^Z(\d+)-(\d+)\.TXT$','tokens','once');
if isempty(match),freqs=cfg.rollHz;else
 hz=str2double(match{2});freqs=unique([cfg.rollHz-hz,cfg.rollHz,cfg.rollHz+hz]);
end
assert(all(freqs>0 & freqs<cfg.fs/2),'无效边带');
end
function path=locateOne(root,filename)
files=[dir(fullfile(root,'**','*.TXT'));dir(fullfile(root,'**','*.txt'));...
 dir(fullfile(root,'*.TXT'));dir(fullfile(root,'*.txt'))];
match=files(arrayfun(@(x)strcmpi(x.name,filename),files));
if ~isempty(match)
 paths=string(fullfile({match.folder},{match.name}));[~,ii]=unique(lower(paths),'stable');match=match(ii);
end
assert(numel(match)==1,'应当找到一个%s，实际%d',filename,numel(match));
path=fullfile(match.folder,match.name);
end
function T=addrow(T,r)
if isempty(T),T=r;else,T=[T;r];end
end
