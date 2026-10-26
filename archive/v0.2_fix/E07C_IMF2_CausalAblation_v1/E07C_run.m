function E07C_run()
%E07C_RUN 只用真实Z3-1数据，检验IMF2硬分类翻转的实际影响。
% 目标：每个alpha分解一次，然后固定同一组modes，仅切换IMF2的重构决策。
% 既不重新运行FVIM优化，也不访问Z3-2或姿态真值。
clc; close all;
cfg=E07B_config();
root=fileparts(mfilename('fullpath'));
outputDir=fullfile(root,'results_E07C');
if ~isfolder(outputDir), mkdir(outputDir); end
assert(isfile(cfg.inputTrain),'找不到Z3-1训练MAT：%s',cfg.inputTrain);
d=load(cfg.inputTrain,'X','baseline','t');
assert(all(isfield(d,{'X','baseline','t'})),'输入MAT应包含X、baseline、t');
X=double(d.X); t=double(d.t(:));
assert(size(X,2)==2 && numel(t)==size(X,1) && all(isfinite(X(:))),...
    '双通道输入形状、时间或有限性不正确');

% 精确复用E07-B先前训练中出现的两组alpha，无新的参数优化。
alphaList=[292.9731,306.1287];
expectedOriginalClass=[1,0];  % 1混叠小波、0丢弃，来自先前CSV核对
K=8;
records=table(); imfRecords=table(); solverRecords=table();
sideF=cfg.trainSidebandHz;
freqs=cfg.trainComponentsHz;
base=E07B_candidate_metrics(X,d.baseline,t,freqs,sideF,cfg);

for q=1:numel(alphaList)
    alpha=alphaList(q);
    fprintf('\n[E07-C %d/2] 对双通道联合MVMD分解一次: K=%d, alpha=%.4f\n',q,K,alpha);
    tic;
    [modes,centers,info]=gyro_mvmd(X,cfg.fs,K,alpha,...
       cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
    elapsed=toc;
    solverRecords=[solverRecords;table(alpha,logical(info.converged),...
       info.iterations,info.AC_reconRelative,elapsed,...
       'VariableNames',{'alpha','converged','iterations','acReconRelative','elapsedSeconds'})]; %#ok<AGROW>
    writetable(solverRecords,fullfile(outputDir,'E07C_solver.csv'));
    if ~info.converged || info.AC_reconRelative>cfg.joint.maxAcReconstructionLoss
        error('E07C:invalidDecomposition',...
            'alpha=%.4f 未收敛或AC重构残差超标。不能据此分析分类切换。',alpha);
    end

    % 保留原始分类策略及所有非IMF2重构处理。
    [Ynatural,details]=gyro_reconstruct(modes,centers,X,cfg,'paper');
    class2=details.class_0drop_1mix_2keep(2);
    if class2~=expectedOriginalClass(q)
        error('E07C:classMismatch',...
            ['alpha=%.4f的IMF2实际分类=%d，预期=%d。请核对MAT、MVMD版本、' ...
             'SampEn版本和config，不得强行解释。'],...
            alpha,class2,expectedOriginalClass(q));
    end
    details.alpha=repmat(alpha,height(details),1);
    imfRecords=[imfRecords;details]; %#ok<AGROW>
    writetable(imfRecords,fullfile(outputDir,'E07C_all_IMF.csv'));

    % 重构仅修改IMF2：两通道各按工程原Haar小波函数生成混叠分量。
    % 不修改IMF1和IMF3～8、不改变已求得的modes、中心频率、均值、其它状态。
    mixed2=zeros(size(X));
    for ch=1:2
        mixed2(:,ch)=gyro_haar(modes(:,ch,2),...
            cfg.wavelet.levels,cfg.wavelet.thresholdScale);
    end
    if class2==1
        Ymixed=Ynatural;
        Ydrop=Ynatural-mixed2;
    elseif class2==0
        Ydrop=Ynatural;
        Ymixed=Ynatural+mixed2;
    else
        error('E07C:unexpectedClass','IMF2本次应仅为mixed/drop。');
    end
    % 检查强制差分恰为IMF2混叠小波贡献，避免额外改变分解本身。
    assert(norm((Ymixed-Ydrop)-mixed2,'fro')/max(norm(mixed2,'fro'),eps)<1e-9,...
       'IMF2切换一致性自检失败');

    for policy=1:2
        if policy==1, Y=Ymixed; name="IMF2_mixed_wavelet"; else
            Y=Ydrop; name="IMF2_dropped"; end
        m=E07B_candidate_metrics(X,Y,t,freqs,sideF,cfg);
        gain=m.nonTargetSuppressionDb(:)-base.nonTargetSuppressionDb(:);
        row=table(string('Z3-1'),K,alpha,string(name),class2,...
            m.nonTargetSuppressionDb(1),m.nonTargetSuppressionDb(2),...
            gain(1),gain(2),...
            m.minSidebandRetentionPct,m.maxSidebandRetentionPct,...
            m.maxAmplitudeDeviationPct,m.maxPhaseShiftDeg,...
            m.maxRelativePhaseShiftDeg,logical(m.sidebandPass),...
            'VariableNames',{'caseName','K','alpha','IMF2Policy','originalIMF2Class',...
            'nonTargetDbA','nonTargetDbS','extraVsFrozenDbA','extraVsFrozenDbS',...
            'minRetentionPct','maxRetentionPct','maxAmplitudeDeviationPct',...
            'maxPhaseShiftDeg','maxRelativePhaseDeg','sidebandPass'});
        records=[records;row]; %#ok<AGROW>
    end
    writetable(records,fullfile(outputDir,'E07C_cross_reconstruction.csv'));
end

% 固定分解、只改分类的因果消融；并检查固定分类时改变alpha造成的差异。
getrow=@(a,label) records(abs(records.alpha-a)<1e-5 & ...
    records.IMF2Policy==string(label),:);
a=alphaList(1);b=alphaList(2);
ma=getrow(a,'IMF2_mixed_wavelet'); da=getrow(a,'IMF2_dropped');
mb=getrow(b,'IMF2_mixed_wavelet'); db=getrow(b,'IMF2_dropped');
assert(height(ma)==1 && height(da)==1 && height(mb)==1 && height(db)==1);
impact=table(...
    [a;b], ...
    [da.nonTargetDbA-ma.nonTargetDbA;db.nonTargetDbA-mb.nonTargetDbA],...
    [da.nonTargetDbS-ma.nonTargetDbS;db.nonTargetDbS-mb.nonTargetDbS],...
    [da.sidebandPass;db.sidebandPass],...
    'VariableNames',{'alpha','dropInsteadOfMixedGainDbA','dropInsteadOfMixedGainDbS','dropSidebandPass'});
writetable(impact,fullfile(outputDir,'E07C_IMF2_causal_delta.csv'));

% 如果使用完全相同的重构规则，alpha改变是否仍导致显著差别？
matched=table(string({'IMF2_mixed_wavelet';'IMF2_dropped'}),...
    [mb.nonTargetDbA-ma.nonTargetDbA;db.nonTargetDbA-da.nonTargetDbA],...
    [mb.nonTargetDbS-ma.nonTargetDbS;db.nonTargetDbS-da.nonTargetDbS],...
    'VariableNames',{'sameIMF2Policy','alpha306Minus293DbA','alpha306Minus293DbS'});
writetable(matched,fullfile(outputDir,'E07C_same_policy_alpha_delta.csv'));

% 这是原因诊断而非新滤波创新通过。固定规则事前写定：两种alpha内部
% drop较mixed在两路均改善>=1dB，且drop均保留边带，就认为硬分类影响值得优化。
classificationImpact=all(impact.dropInsteadOfMixedGainDbA>=1) && ...
    all(impact.dropInsteadOfMixedGainDbS>=1) && all(impact.dropSidebandPass);
smallAlphaAfterMatched=all(abs(matched.alpha306Minus293DbA)<=0.5) && ...
    all(abs(matched.alpha306Minus293DbS)<=0.5);
fid=fopen(fullfile(outputDir,'E07C_decision.txt'),'w','n','UTF-8');
assert(fid>0,'无法创建decision文件');
fprintf(fid,'E07-C fixed-decomposition IMF2 cross-reconstruction (Z3-1 TRAIN ONLY)\n');
fprintf(fid,'Two MVMD solves; no FVIM search; no Z3-2 data accessed.\n');
fprintf(fid,'IMF2 hard classification is different: alpha293=mixed; alpha306=dropped.\n');
fprintf(fid,'Same-decomposition IMF2 drop-over-mixed change in dB:\n');
for p=1:height(impact)
    fprintf(fid,'  alpha %.4f: channel a %.6f dB, channel s %.6f dB; sidebandPass=%d\n',...
        impact.alpha(p),impact.dropInsteadOfMixedGainDbA(p),...
        impact.dropInsteadOfMixedGainDbS(p),impact.dropSidebandPass(p));
end
fprintf(fid,'Classification effect passes diagnostic criterion: %d\n',classificationImpact);
fprintf(fid,'Alpha-only difference with a fixed class <=0.5dB in both channels: %d\n',smallAlphaAfterMatched);
if classificationImpact && smallAlphaAfterMatched
    fprintf(fid,'INTERPRETATION: strong evidence the discontinuous IMF2 decision explains most noise-proxy gap.\n');
    fprintf(fid,'NEXT: develop a continuous paired-mode weighting rule ON TRAIN DATA; evaluate on NEW separate records.\n');
elseif classificationImpact
    fprintf(fid,'INTERPRETATION: IMF2 class switch matters, but decomposition-alpha difference also matters.\n');
    fprintf(fid,'NEXT: quantify joint IMF classification stability; avoid asserting class switch explains all differences.\n');
else
    fprintf(fid,'INTERPRETATION: class-switch explanation not fully supported; recheck numeric decomposition/metrics.\n');
end
fprintf(fid,'Not a true clean-signal SNR, body Omega_x/Omega_y or attitude RMSE.\n');
fclose(fid);
fprintf('\nE07-C完成。请返回E07C_cross_reconstruction.csv、E07C_IMF2_causal_delta.csv、E07C_same_policy_alpha_delta.csv、E07C_decision.txt。\n');
end
