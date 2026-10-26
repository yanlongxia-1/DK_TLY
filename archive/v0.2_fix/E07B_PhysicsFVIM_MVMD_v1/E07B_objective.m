function [entropyFitness,jointFitness,entropyDetail,metrics] = E07B_objective(X,Y,modes,t,frequencies,sidebands,cfg)
% entropy-only: 仿照早期论文工程代理，模态最小包络熵；非论文可执行适应度的确切复刻。
% joint: 固定的多项物理保真+非目标残差代理；优选为最小化。
K=size(modes,3);
ent=nan(K,2);
for k=1:K
    for c=1:2
        env=abs(gyro_analytic(modes(:,c,k)));
        p=env/(sum(env)+eps);
        ent(k,c)=-sum(p.*log(p+eps))/log(numel(p));
    end
end
entropyFitness=mean(min(ent,[],1));
energy=squeeze(sum(sum(modes.^2,1),2));
energyShare=energy/(sum(energy)+eps);
if any(energyShare<1e-5),entropyFitness=entropyFitness+1;end
metrics=E07B_candidate_metrics(X,Y,t,frequencies,sidebands,cfg);
w=cfg.joint.weights;
noiseTerm=mean(metrics.nonTargetResidualRatio.^2);
ampTerm=mean((metrics.retentionPct(:)-100).^2)/cfg.joint.amplitudeScalePct^2;
ph=metrics.phaseAbsDeg(metrics.phaseValid);
phaseTerm=mean((ph/cfg.joint.phaseScaleDeg).^2);
relativeTerm=mean((metrics.relativePhaseChangeDeg/cfg.joint.phaseScaleDeg).^2);
jointFitness=w(1)*noiseTerm+w(2)*ampTerm+w(3)*phaseTerm+w(4)*relativeTerm;
entropyDetail=struct('noiseTerm',noiseTerm,'amplitudeTerm',ampTerm,...
    'phaseTerm',phaseTerm,'relativeTerm',relativeTerm,'entropyMeanMin',entropyFitness);
end
