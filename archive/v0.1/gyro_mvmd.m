function [modes, omegaHz, info] = gyro_mvmd(X, fs, K, alpha, maxIter, tol, tau)
% Self-contained two/multichannel VMD with shared mode center frequencies.
% Inputs: N x C, modes: N x C x K. Offline/mirrored extension.
% FFT uses frequency in cycles/sample; alpha here follows THIS normalization.
X=double(X); [N,C]=size(X);
if N<100 || K<1 || alpha<=0, error('Invalid MVMD input/parameter'); end
half=floor(N/2);
Xm=[flipud(X(1:half,:)); X; flipud(X(end-half+1:end,:))];
if mod(size(Xm,1),2)~=0, Xm=[Xm; Xm(end,:)]; end
M=size(Xm,1);
means=mean(Xm,1);
Xm=Xm-means;
F=fftshift(fft(Xm,[],1),1);
F(1:M/2,:)=0; % one-sided representation, no analytic doubling
f=((0:M-1)'/M)-0.5;
U=complex(zeros(M,C,K));
omega=linspace(0,0.5,K+1); omega=omega(1:K);
if K==1, omega=0.005; end
lambda=complex(zeros(M,C));
converged=false;
for it=1:maxIter
    Uold=U;
    for k=1:K
        other=sum(U,3)-U(:,:,k);
        denom=1+2*alpha*(f-omega(k)).^2;
        U(:,:,k)=(F-other-lambda/2)./denom;
        pow=sum(abs(U(:,:,k)).^2,2);
        d=sum(pow(M/2+1:end));
        if d>eps
            omega(k)=sum(f(M/2+1:end).*pow(M/2+1:end))/d;
        end
        omega(k)=min(0.5,max(0,omega(k)));
    end
    lambda=lambda+tau*(sum(U,3)-F);
    delta=norm(U(:)-Uold(:))/(norm(Uold(:))+eps);
    if delta<tol && it>10
        converged=true; break
    end
end
omegaHz=omega(:)*fs;
[omegaHz,order]=sort(omegaHz);
U=U(:,:,order);
modes=zeros(N,C,K);
for k=1:K
    % Rebuild Hermitian Fourier spectrum from shared positive-frequency modes.
    spec=zeros(M,C);
    spec(M/2+1:end,:)=U(M/2+1:end,:,k);
    spec(2:M/2,:)=conj(flipud(U(M/2+2:end,:,k)));
    fullWave=real(ifft(ifftshift(spec,1),[],1));
    modes(:,:,k)=fullWave(half+1:half+N,:);
end
% Measured DC is carried separately by the caller; don't assign it to a noise IMF.
info.iterations=it;
info.converged=converged;
info.delta=delta;
info.dc=mean(X,1);
info.reconRelative=norm(reshape(sum(modes,3),[],1)-reshape(X-mean(X,1),[],1))/max(norm(X(:)),eps);
end
