function [s,diag]=E07D_score(X,Y,modes,t,freqs,sideF,cfg)
%E07D_SCORE 同一原始信号下固定物理代理指标；训练选择不碰外部留出数据。
m=E07B_candidate_metrics(X,Y,t,freqs,sideF,cfg);
w=cfg.joint.weights;
noise=mean(m.nonTargetResidualRatio(:).^2);
amp=mean(((m.retentionPct(:)-100)/cfg.joint.amplitudeScalePct).^2);
valid=m.phaseAbsDeg(m.phaseValid);
if isempty(valid),valid=0;end
phase=mean((valid/cfg.joint.phaseScaleDeg).^2);
rel=mean((m.relativePhaseChangeDeg(:)/cfg.joint.phaseScaleDeg).^2);
parts=[noise amp phase rel];
% 不满足事前保真窗的候选，不得仅凭过度降噪取得适应度优势。
ret=m.retentionPct(:);
low=max(0,cfg.pass.minRetentionPct-ret);
high=max(0,ret-cfg.pass.maxRetentionPct);
validPhase=m.phaseAbsDeg(m.phaseValid);
phaseViolation=max(0,validPhase-cfg.pass.maxSidebandPhaseShiftDeg);
constraintPenalty=cfg.joint.constraintPenalty*(sum(low.^2+high.^2)+...
    sum((phaseViolation/cfg.joint.phaseScaleDeg).^2));
s=sum(w.*parts)+constraintPenalty;
if ~isfinite(s),s=Inf;end
diag.metrics=m;
diag.parts=parts;
diag.constraintPenalty=constraintPenalty;
end
