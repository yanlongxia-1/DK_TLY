function [summary,details] = E08J_physics_metrics(X,Y,t,freqs,cfg)
% E08-J: stable relative fidelity for strong spectral components and bounded
% complex-coefficient distortion for weak components, jointly across a/s.
% This is a physics-informed PROXY relative to raw measurements, not truth.
assert(isequal(size(X),size(Y)) && size(X,2)==2, 'Require two equal-length channels');
t=t(:);N=size(X,1);assert(numel(t)==N && N>numel(freqs)*3);
freqs=unique(freqs(:)');assert(all(freqs>0 & freqs<cfg.fs/2));
F=ones(N,1);
for k=1:numel(freqs)
    F=[F cos(2*pi*freqs(k)*t) sin(2*pi*freqs(k)*t)]; %#ok<AGROW>
end
bX=F\X;bY=F\Y;
R= X-F*bX; S=Y-F*bY;
resX=sqrt(mean(R.^2,1));resY=sqrt(mean(S.^2,1));
proxyDb=20*log10(max(resX,eps)./max(resY,eps));
M=numel(freqs); Cx=complex(zeros(2,M));Cy=Cx;
for j=1:M
    Cx(:,j)=complex(bX(2*j,:),bX(2*j+1,:)).';
    Cy(:,j)=complex(bY(2*j,:),bY(2*j+1,:)).';
end
AmpX=abs(Cx);AmpY=abs(Cy);
ampReference=max(AmpX,[],2);assert(all(ampReference>1e-12), 'Cannot determine nonzero raw reference amplitude');
strong=AmpX>=cfg.policy.strongFraction * repmat(ampReference,1,M);
retention=100*AmpY./max(AmpX,eps);
phaseChange=abs(wrapD(angle(Cy./Cx)*180/pi));
complexError=abs(Cy-Cx);
normalizedComplexError=complexError./repmat(ampReference,1,M);
strongPass=abs(retention-100)<=cfg.policy.strongAmpTolerancePct & ...
           phaseChange<=cfg.policy.strongPhaseToleranceDeg;
weakPass=normalizedComplexError<=cfg.policy.weakComplexToleranceFraction;
cellPass= (strong&strongPass)|(~strong&weakPass);
relativePhase=zeros(1,M);relativePass=true(1,M);oldRelativePass=true(1,M);
for j=1:M
  raw=angle(Cx(2,j))-angle(Cx(1,j));
  out=angle(Cy(2,j))-angle(Cy(1,j));
  relativePhase(j)=abs(wrapD((out-raw)*180/pi));
  % Only meaningful when BOTH raw coefficients are strong; weak components
  % are already checked by absolute complex error in EACH channel.
  if all(strong(:,j))
      relativePass(j)=relativePhase(j)<=cfg.policy.strongRelativePhaseToleranceDeg;
  end
  oldRelativePass(j)=relativePhase(j)<=cfg.pass.maxSidebandPhaseShiftDeg;
end
% Historical E08H/E08I logic: 1% amplitude everywhere + phase on features
% with amp >= 1% of the global maximum + relative phase for all frequencies.
oldPhaseValid=AmpX>=.01*max(AmpX(:));
oldCellPass=retention>=cfg.pass.minRetentionPct & retention<=cfg.pass.maxRetentionPct & ...
  ((~oldPhaseValid)|(phaseChange<=cfg.pass.maxSidebandPhaseShiftDeg));
summary.newPass=all(cellPass(:))&&all(relativePass);
summary.oldPass=all(oldCellPass(:))&&all(oldRelativePass);
summary.minimumStrongRetention=NaN;summary.maximumStrongRetention=NaN;
if any(strong(:))
  summary.minimumStrongRetention=min(retention(strong));
  summary.maximumStrongRetention=max(retention(strong));
end
summary.maximumWeakComplexFraction=0;
if any(~strong(:))
  summary.maximumWeakComplexFraction=max(normalizedComplexError(~strong));
end
summary.minimumProxyDb=min(proxyDb);
summary.meanProxyDb=mean(proxyDb);
summary.maximumStrongPhaseDeg=0;
if any(strong(:)),summary.maximumStrongPhaseDeg=max(phaseChange(strong));end
summary.maximumStrongRelativePhaseDeg=0;
if any(arrayfun(@(j)all(strong(:,j)),1:M))
  mask=arrayfun(@(j)all(strong(:,j)),1:M);
  summary.maximumStrongRelativePhaseDeg=max(relativePhase(mask));
end
channels=repelem((1:2)',M);fcol=repmat(freqs(:),2,1);
strongColumn=reshape(strong.',[],1);
% Column ordering matches channels: [ch1 f1..fM, ch2 f1..fM].
idx=[1:M,M+1:2*M]; %#ok<NASGU>
rawColumn=reshape(AmpX.',[],1);
refColumn=repelem(ampReference,M);
retColumn=reshape(retention.',[],1);
phaseColumn=reshape(phaseChange.',[],1);
relColumn=repmat(relativePhase(:),2,1);
errColumn=reshape(complexError.',[],1);
normErrColumn=reshape(normalizedComplexError.',[],1);
passColumn=reshape(cellPass.',[],1);
relPassColumn=repmat(relativePass(:),2,1);
proxyColumn=repelem(proxyDb(:),M);
details=table(channels,fcol,rawColumn,refColumn,logical(strongColumn),...
 retColumn,phaseColumn,relColumn,errColumn,normErrColumn,logical(passColumn),...
 logical(relPassColumn),proxyColumn, 'VariableNames',{'channel','frequencyHz',...
 'rawAmplitude','strongReferenceAmplitude','isStrong','retentionPct',...
 'phaseChangeDeg','relativePhaseChangeDeg','complexCoefficientError',...
 'normalizedComplexError','newCellPass','newRelativePhasePass','proxyDb'});
end
function w=wrapD(x),w=mod(x+180,360)-180;end
