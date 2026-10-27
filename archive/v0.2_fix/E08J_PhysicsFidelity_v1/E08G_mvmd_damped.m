function [modes, omegaHz, info, snaps] = E08G_mvmd_damped(X,fs,K,alpha,maxIter,tol,tau,rollHz,eta,checkpoints)
% E08-G / E08F solver-equation reproduction; ONLY change: center update relaxation.
% eta=1 reproduces original update; 0<eta<1 dampens center movement.
% Changes in eta do NOT change the fixed-point center-frequency equation.
% This solver and its output are diagnostic until separately validated.
assert(isfinite(eta)&&eta>0&&eta<=1,'Require 0<eta<=1');
X=double(X);[N,C]=size(X);
assert(N>=100&&C==2&&K>=1&&K==round(K)&&alpha>0&&fs>0);
assert(all(isfinite(X(:)))&&maxIter>=1&&tol>0&&tau>=0&&rollHz>0&&rollHz<fs/2);
h=floor(N/2);Xm=[flipud(X(1:h,:));X;flipud(X(end-h+1:end,:))];
if mod(size(Xm,1),2),Xm=[Xm;Xm(end,:)];end
M=size(Xm,1);half=M/2;Xm=Xm-mean(Xm,1);
F=fftshift(fft(Xm,[],1),1);F(1:half,:)=0;
f=((0:M-1)'/M)-.5;pos=(half+1):M;
U=complex(zeros(M,C,K));SumU=complex(zeros(M,C));lambda=complex(zeros(M,C));
if K==1,seedHz=rollHz;
elseif K==2,seedHz=[rollHz min(fs*.44,max(2*rollHz,rollHz+3))];
else
 low2=min(fs*.40,max(2*rollHz,rollHz+3));
 highStart=max(40,low2+4);
 seedHz=[rollHz low2 linspace(highStart,fs*.46,K-2)];
end
omega=sort(seedHz(:)')/fs;
snaps=struct('iteration',{},'modes',{},'centerHz',{},'delta',{},'centerDeltaHz',{},'centerStepHz',{},'modeLocalDelta',{},'modeGlobalDeltaShare',{},'modeEnergyShare',{},'ACreconRelative',{});
countStable=0;converged=false;deltaHistory=zeros(maxIter,1);centerHistory=zeros(maxIter,1);
for it=1:maxIter
 Uold=U;oOld=omega;centerGaps=zeros(K,1);
 for k=1:K
  SumU=SumU-U(:,:,k);
  denom=1+2*alpha*(f-omega(k)).^2;
  Uk=(F-SumU-lambda/2)./denom;Uk(1:half,:)=0;
  U(:,:,k)=Uk;SumU=SumU+Uk;
  power=sum(abs(Uk(pos,:)).^2,2);ptotal=sum(power);
  if ptotal>eps
   target=sum(f(pos).*power)/ptotal;
   target=min(.5,max(0,target));
   centerGaps(k)=abs(target-omega(k))*fs;
   omega(k)=(1-eta)*omega(k)+eta*target;
  end
  omega(k)=min(.5,max(0,omega(k)));
 end
 lambda=lambda+tau*(SumU-F);
 D=U-Uold;
 delta=norm(D(:))/max(norm(Uold(:)),eps);
 centerStep=max(abs(omega-oOld))*fs;
 % Compare UNRELAXED target-frequency gap: eta<1 does not get a looser center test.
 centerDelta=max(centerGaps);
 deltaHistory(it)=delta;centerHistory(it)=centerDelta;
 if ~isfinite(delta)||~isfinite(centerDelta),error('Nonfinite solver state');end
 if delta<tol && centerDelta<max(fs*tol,.001) && it>=12
  countStable=countStable+1;
 else
  countStable=0;
 end
 if countStable>=3,converged=true;end
 if any(checkpoints==it)||converged||it==maxIter
  ener=zeros(K,1);local=zeros(K,1);change=zeros(K,1);
  for j=1:K
   old=Uold(:,:,j);u=U(:,:,j);d=u-old;
   ener(j)=sum(abs(u(:)).^2);
   local(j)=norm(d(:))/max(norm(old(:)),eps);
   change(j)=sum(abs(d(:)).^2);
  end
  [centSorted,order]=sort(omega(:)*fs);
  [mm,recon]=unmirror(U(:,:,order),X,h,half,M,N,C,K,pos);
  sn.iteration=it;sn.modes=mm;sn.centerHz=centSorted;
  sn.delta=delta;sn.centerDeltaHz=centerDelta;sn.centerStepHz=centerStep;
  sn.modeLocalDelta=local(order);
  sn.modeGlobalDeltaShare=change(order)/(sum(change)+eps);
  sn.modeEnergyShare=ener(order)/(sum(ener)+eps);
  sn.ACreconRelative=recon;
  snaps(end+1)=sn; %#ok<AGROW>
 end
 if converged,break;end
end
[modes,omegaHz]=deal(snaps(end).modes,snaps(end).centerHz);
info.iterations=it;info.converged=converged;info.delta=delta;
info.centerDeltaHz=centerDelta;info.centerStepHz=centerStep;
info.ACreconRelative=snaps(end).ACreconRelative;
info.deltaHistory=deltaHistory(1:it);info.centerHistoryHz=centerHistory(1:it);
info.eta=eta;
end
function [modes,acerr]=unmirror(U,X,h,half,M,N,C,K,pos)
modes=zeros(N,C,K);
for k=1:K
 spec=zeros(M,C);spec(pos,:)=U(pos,:,k);
 spec(2:half,:)=conj(flipud(U(half+2:end,:,k)));
 wave=real(ifft(ifftshift(spec,1),[],1));
 modes(:,:,k)=wave(h+1:h+N,:);
end
AC=X-mean(X,1);
R=AC-sum(modes,3);
acerr=norm(R(:))/max(norm(AC(:)),eps);
end
