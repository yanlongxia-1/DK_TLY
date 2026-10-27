function [modes, omegaHz, info] = gyro_mvmd_conditioned(X, fs, K, alpha, maxIter, tol, tau, rollHz, useConditioned)
% Joint-channel, mirrored MVMD (noise-tolerant with tau=0).
% Adapted from the E08A baseline: ONLY initialization is modified.
% This is not claimed as the original paper author FVIM implementation.
% X: N-by-C, returned modes: N-by-C-by-K. Same center Hz across channels.
% IMPORTANT: tau=0 solves a noise-tolerant relaxed reconstruction, therefore
% sum(modes) need not equal X-mean(X). Track AC reconstruction loss separately.
% With tau>0 the dual-ascent solver can oscillate; it is intended for diagnosis,
% and NONCONVERGED solutions must NOT be used in final comparisons.
X = double(X); [N,C] = size(X);
if N<100 || C<1 || K<1 || K~=round(K) || alpha<=0 || fs<=0 || ...
        maxIter<1 || tol<=0 || tau<0 || any(~isfinite(X(:)))
    error('gyro_mvmd:badInput','Invalid MVMD inputs or parameters.');
end
h = floor(N/2);
Xm = [flipud(X(1:h,:)); X; flipud(X(end-h+1:end,:))];
if mod(size(Xm,1),2), Xm=[Xm;Xm(end,:)];end
M = size(Xm,1); half=M/2;
Xm = Xm-mean(Xm,1);
F = fftshift(fft(Xm,[],1),1);
F(1:half,:) = 0;
f = ((0:M-1)'/M)-0.5; pos=(half+1):M;
U = complex(zeros(M,C,K));
SumU = complex(zeros(M,C));
% Frequency-informed deterministic INITIALIZATION only; FVIM still optimizes K/alpha.
% The reconstruction equations, tau, convergence tests and limits are unchanged.
if useConditioned
    assert(isfinite(rollHz) && rollHz>0 && rollHz<fs/2);
    if K==1
        seedHz=rollHz;
    elseif K==2
        seedHz=[rollHz min(fs*0.44,max(2*rollHz,rollHz+3))];
    else
        low2=min(fs*0.40,max(2*rollHz,rollHz+3));
        highStart=max(40,low2+4);
        seedHz=[rollHz, low2, linspace(highStart,fs*0.46,K-2)];
    end
    omega=sort(seedHz(:)')/fs;
else
    omega=linspace(0,0.5,K+1);omega=omega(1:K);
    if K==1,omega=0.005;end
    seedHz=omega*fs;
end
lambda = complex(zeros(M,C));
converged=false; countStable=0;
deltaHist=zeros(maxIter,1); centerHist=zeros(maxIter,1); residualHist=zeros(maxIter,1);
for it=1:maxIter
    Uold=U; oOld=omega;
    for k=1:K
        % Remove previous k-th mode from running sum, then Gauss-Seidel update.
        SumU=SumU-U(:,:,k);
        denom=1+2*alpha*(f-omega(k)).^2;
        Uk=(F-SumU-lambda/2)./denom;
        Uk(1:half,:)=0;
        U(:,:,k)=Uk;
        SumU=SumU+Uk;
        spectralPower=sum(abs(Uk(pos,:)).^2,2);
        ptotal=sum(spectralPower);
        if ptotal>eps
            omega(k)=sum(f(pos).*spectralPower)/ptotal;
        end
        omega(k)=min(0.5,max(0,omega(k)));
    end
    lambda = lambda+tau*(SumU-F);
    delta=norm(U(:)-Uold(:))/max(norm(Uold(:)),eps);
    odelta=max(abs(omega-oOld))*fs;
    fres=norm(SumU(:)-F(:))/max(norm(F(:)),eps);
    deltaHist(it)=delta;centerHist(it)=odelta;residualHist(it)=fres;
    if ~isfinite(delta) || ~isfinite(odelta) || ~isfinite(fres)
        error('gyro_mvmd:numericalFailure','Nonfinite MVMD iteration state.');
    end
    % Require BOTH spectral-mode and center-frequency stabilization.
    % Three consecutive checks avoid a spurious isolated threshold crossing.
    if delta<tol && odelta<max(fs*tol,0.001) && it>=12
        countStable=countStable+1;
    else
        countStable=0;
    end
    if countStable>=3,converged=true;break;end
end
omegaHz=omega(:)*fs;
[omegaHz,order]=sort(omegaHz);
U=U(:,:,order);
modes=zeros(N,C,K);
for k=1:K
    spec=zeros(M,C);
    spec(pos,:)=U(pos,:,k);
    % Undo fftshift; explicitly enforce Hermitian symmetry.
    spec(2:half,:)=conj(flipud(U(half+2:end,:,k)));
    fullWave=real(ifft(ifftshift(spec,1),[],1));
    modes(:,:,k)=fullWave(h+1:h+N,:);
end
AC=X-mean(X,1);
R=AC-sum(modes,3);
info.iterations=it;
info.converged=converged;
info.delta=deltaHist(it);
info.centerDeltaHz=centerHist(it);
info.spectralResidual=residualHist(it);
info.AC_reconRelative=norm(R(:))/max(norm(AC(:)),eps);
info.reconRelative=info.AC_reconRelative; % compatibility, NEW definition
info.dc=mean(X,1);
info.tau=tau;
info.initialCenterHz=seedHz;
info.conditionedInitialization=logical(useConditioned);
info.status='nonconverged';
if converged,info.status='converged';end
info.deltaHistory=deltaHist(1:it);
info.centerHistoryHz=centerHist(1:it);
info.residualHistory=residualHist(1:it);
end
