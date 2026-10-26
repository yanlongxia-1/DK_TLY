%% E07-A: TRAIN-only small-budget MVMD alpha precheck for later joint FVIM fitness.
% Objective: can changing MVMD alpha improve non-target residual proxy by >=1dB
% in BOTH channels WITHOUT distorting known motion-related sidebands?
% WARNING: not FVIM optimizer; not a proven algorithm innovation or true SNR.
clear; clc; close all;
cfg = E07A_config();
assert(isfile(cfg.inputTrain) && isfile(cfg.inputTest), 'Missing E06B MAT files under inputs/');
if ~isfolder(cfg.resultDir), mkdir(cfg.resultDir); end
tr = load(cfg.inputTrain,'X','baseline','t');
assert(all(isfield(tr,{'X','baseline','t'})),'Training MAT lacks X/baseline/t');
assert(isequal(size(tr.X),size(tr.baseline)) && size(tr.X,2)==2,'Data dimension mismatch');
trainF=[2,3,4]; sideF=[2,4];
trainTable = table(); solverTable=table();
% Original frozen alpha baseline is already computed, avoiding any refitting.
baseSummary = calcSummary(tr.X,tr.baseline,tr.baseline,tr.t,trainF,sideF,cfg,'Z3-1',cfg.mvmd.alphaBaseline,'baseline');
trainTable = [trainTable;baseSummary];
selectedAlpha=NaN; selectedScore=-Inf;
for a=cfg.mvmd.alphaCandidates
    fprintf('[E07-A TRAIN] alpha=%.8g, K=%d; one joint MVMD\n',a,cfg.mvmd.K);
    tic;
    [modes,centerHz,info] = gyro_mvmd(tr.X,cfg.fs,cfg.mvmd.K,a, ...
        cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
    seconds=toc;
    solverTable=[solverTable;table(string('Z3-1'),a,logical(info.converged), ...
        info.iterations,info.AC_reconRelative,seconds, ...
        'VariableNames',{'caseName','alpha','converged','iterations','acReconLoss','elapsedSeconds'})];
    if ~info.converged
        fprintf('  Nonconverged: skip. No extension or retuning.\n');
        continue;
    end
    [candidate,~]=gyro_reconstruct(modes,centerHz,tr.X,cfg,'paper');
    s=calcSummary(tr.X,tr.baseline,candidate,tr.t,trainF,sideF,cfg,'Z3-1',a,'candidate');
    trainTable=[trainTable;s];
    pass = all(s.sidebandPass) && all(s.extraSuppressionDb>=cfg.pass.minExtraSuppressionDbBothChannels);
    score=min(s.extraSuppressionDb);
    fprintf('  additional suppression a/s=%.3f/%.3f dB, sidebandPass=%d, pass=%d\n', ...
        s.extraSuppressionDb(1),s.extraSuppressionDb(2),all(s.sidebandPass),pass);
    if pass && score>selectedScore
        selectedScore=score;
        selectedAlpha=a;
    end
end
writetable(trainTable,fullfile(cfg.resultDir,'E07A_train_objective.csv'));
writetable(solverTable,fullfile(cfg.resultDir,'E07A_solver.csv'));
if ~isfinite(selectedAlpha)
    decision='STOP: no fixed-alpha candidate passed train. Do not test/tune Z3-2. Parameter-only change not supported.';
    fid=fopen(fullfile(cfg.resultDir,'E07A_decision.txt'),'w');
    fprintf(fid,'%s\n',decision);
    fprintf(fid,'MVMD baseline remains as-is. Proceed to a substantively different, explicitly specified FVIM objective or reconstruction model.\n');
    fclose(fid);
    disp(decision); return;
end
% One HELD-OUT test (never used in selecting alpha).
te=load(cfg.inputTest,'X','baseline','t');
assert(all(isfield(te,{'X','baseline','t'})),'Test MAT lacks X/baseline/t');
fTest=[1,3,5]; sideTest=[1,5];
tic;
[modes,centers,info]=gyro_mvmd(te.X,cfg.fs,cfg.mvmd.K,selectedAlpha, ...
    cfg.mvmd.maxIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
seconds=toc;
solverTable=[solverTable;table(string('Z3-2'),selectedAlpha,logical(info.converged), ...
    info.iterations,info.AC_reconRelative,seconds, ...
    'VariableNames',{'caseName','alpha','converged','iterations','acReconLoss','elapsedSeconds'})];
writetable(solverTable,fullfile(cfg.resultDir,'E07A_solver.csv'));
if info.converged
    [testCandidate,~]=gyro_reconstruct(modes,centers,te.X,cfg,'paper');
    testTable=calcSummary(te.X,te.baseline,testCandidate,te.t,fTest,sideTest,cfg,'Z3-2',selectedAlpha,'heldout');
    writetable(testTable,fullfile(cfg.resultDir,'E07A_heldout_objective.csv'));
    testPass=all(testTable.sidebandPass) && all(testTable.extraSuppressionDb>=cfg.pass.minExtraSuppressionDbBothChannels);
else
    testPass=false;
end
if testPass
    decision=sprintf('PROVISIONAL PASS on held-out Z3-2: alpha %.10g. Still NOT FVIM innovation or true SNR proof.',selectedAlpha);
else
    decision=sprintf('STOP after held-out Z3-2: alpha %.10g did not meet constraints. Do not retune on test.',selectedAlpha);
end
fid=fopen(fullfile(cfg.resultDir,'E07A_decision.txt'),'w');
fprintf(fid,'%s\n',decision); fclose(fid);
disp(decision);

function summary=calcSummary(raw,baseline,candidate,t,targetF,sideF,cfg,caseName,alpha,role)
    N=size(raw,1); t=t(:);
    assert(size(raw,2)==2 && numel(t)==N,'N by 2 data/time mismatch');
    sideRet=zeros(2,numel(sideF)); sidePhase=zeros(2,numel(sideF));
    db=zeros(2,1); perPass=false(2,numel(sideF));
    for c=1:2
        [rC,~]=fitHarmonics(raw(:,c),t,targetF);
        [~,bResidual]=fitHarmonics(baseline(:,c),t,targetF);
        [cC,cResidual]=fitHarmonics(candidate(:,c),t,targetF);
        db(c)=20*log10(max(rmsPlain(bResidual),eps)/max(rmsPlain(cResidual),eps));
        % Re-express improvement relative to the SAME frozen baseline.
        % baseline candidate yields exactly 0dB.
        for j=1:numel(sideF)
            p=find(abs(targetF-sideF(j))<1e-12,1);
            rawAmp=abs(rC(p));
            sideRet(c,j)=100*abs(cC(p))/max(rawAmp,eps);
            sidePhase(c,j)=wrap180(localPhase(cC(p))-localPhase(rC(p)));
            perPass(c,j)=sideRet(c,j)>=cfg.pass.minSidebandRetentionPct && ...
                sideRet(c,j)<=cfg.pass.maxSidebandRetentionPct && ...
                abs(sidePhase(c,j))<=cfg.pass.maxSidebandPhaseShiftDeg;
        end
    end
    % Construct typed table columns directly.  MATLAB cell2table may already
    % convert homogeneous scalar cells to double/logical; calling cell2mat
    % on those columns causes 'input must be cell' on some releases.
    summary=table( ...
        repmat(string(caseName),2,1), ...
        repmat(string(role),2,1), ...
        repmat(double(alpha),2,1), ...
        string({'a';'s'}), ...
        db, ...
        min(sideRet,[],2), ...
        max(sideRet,[],2), ...
        max(abs(sidePhase),[],2), ...
        all(perPass,2), ...
        'VariableNames',{'caseName','role','alpha','channel', ...
        'extraSuppressionDb','minSidebandRetentionPct','maxSidebandRetentionPct', ...
        'maxSidebandPhaseShiftAbsDeg','sidebandPass'});
end
function [coef,residual]=fitHarmonics(y,t,f)
    A=ones(numel(t),1);
    for k=1:numel(f), A=[A,cos(2*pi*f(k)*t),sin(2*pi*f(k)*t)];end %#ok<AGROW>
    beta=A\y;
    coef=complex(beta(2:2:end),beta(3:2:end));
    residual=y-A*beta;
end
function a=localPhase(z),a=atan2d(imag(z),real(z));end
function a=wrap180(x),a=mod(x+180,360)-180;end
function x=rmsPlain(y),x=sqrt(mean(y(:).^2));end
