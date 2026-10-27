function m=E07B_candidate_metrics(X,Y,t,freqs,sideF,cfg)
% 只用相同的真实原始测量作频域代理；不依赖转台真值/标定值。
% C+i*S是COS/SIN系数对；并不直接当成可乘除的解析信号。
assert(isequal(size(X),size(Y)) && size(X,2)==2,'X/Y必须为同形N×2');
t=t(:); assert(numel(t)==size(X,1),'t长度不匹配');
F=ones(numel(t),1);
for j=1:numel(freqs),F=[F,cos(2*pi*freqs(j)*t),sin(2*pi*freqs(j)*t)];end %#ok<AGROW>
ampPct=zeros(2,numel(sideF));phaseAbs=zeros(2,numel(sideF));
relPhaseChange=zeros(1,numel(sideF));
coeffRaw=complex(zeros(2,numel(sideF)));coeffOut=coeffRaw;
noiseDb=zeros(2,1);noiseRatio=zeros(2,1);
rawAmplitude=zeros(2,numel(sideF));
for ch=1:2
    betaRaw=F\X(:,ch); betaY=F\Y(:,ch);
    rX=X(:,ch)-F*betaRaw;rY=Y(:,ch)-F*betaY;
    rmsX=sqrt(mean(rX.^2));rmsY=sqrt(mean(rY.^2));
    noiseRatio(ch)=rmsY/max(rmsX,eps);
    noiseDb(ch)=20*log10(max(rmsX,eps)/max(rmsY,eps));
    for jj=1:numel(sideF)
        p=find(abs(freqs-sideF(jj))<1e-8,1);
        assert(~isempty(p),'sideF并非freqs子集');
        cX=complex(betaRaw(2*p),betaRaw(2*p+1));
        cY=complex(betaY(2*p),betaY(2*p+1));
        coeffRaw(ch,jj)=cX;coeffOut(ch,jj)=cY;
        rawAmplitude(ch,jj)=abs(cX);
        ampPct(ch,jj)=100*abs(cY)/max(abs(cX),eps);
        phaseAbs(ch,jj)=abs(wrap180(atan2d(imag(cY),real(cY)) - atan2d(imag(cX),real(cX))));
    end
end
for jj=1:numel(sideF)
    r1=atan2d(imag(coeffRaw(2,jj)),real(coeffRaw(2,jj)))-atan2d(imag(coeffRaw(1,jj)),real(coeffRaw(1,jj)));
    r2=atan2d(imag(coeffOut(2,jj)),real(coeffOut(2,jj)))-atan2d(imag(coeffOut(1,jj)),real(coeffOut(1,jj)));
    relPhaseChange(jj)=abs(wrap180(r2-r1));
end
% 只有幅值显著的频点才对相位进行稳定约束（两路共享可比性）。
phaseValid=rawAmplitude>=cfg.joint.phaseFloorRawAmpFraction*max(rawAmplitude(:));
if ~any(phaseValid(:)), phaseValid(:)=true; end
m.retentionPct=ampPct;
m.phaseAbsDeg=phaseAbs;
m.relativePhaseChangeDeg=relPhaseChange;
m.phaseValid=phaseValid;
m.nonTargetResidualRatio=noiseRatio;
m.nonTargetSuppressionDb=noiseDb;
m.maxAmplitudeDeviationPct=max(abs(ampPct(:)-100));
m.maxPhaseShiftDeg=max(phaseAbs(phaseValid));
m.maxRelativePhaseShiftDeg=max(relPhaseChange);
m.minSidebandRetentionPct=min(ampPct(:));
m.maxSidebandRetentionPct=max(ampPct(:));
m.sidebandPass=all(ampPct(:)>=cfg.pass.minRetentionPct & ampPct(:)<=cfg.pass.maxRetentionPct) && ...
    all(phaseAbs(phaseValid)<=cfg.pass.maxSidebandPhaseShiftDeg);
end
function out=wrap180(x),out=mod(x+180,360)-180;end
