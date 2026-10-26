function yhat = E05_KF_fixed(y,fs,f0,trainIdx,qOscOverR,qBiasOverR)
% Train-only noise estimate; full sequence recursive causal oscillator KF.
% No look-ahead mean on the validation interval.
[N,C] = size(y);
if C ~= 2, error('Expected two channels [a,s].');end
yhat = zeros(N,C);
fit = gyro_sinefit(y(trainIdx,:),fs,f0);
Rtrain = max(var(fit.residual,0,1),1e-12);
ang=2*pi*f0/fs;
A=[cos(ang),-sin(ang),0;sin(ang),cos(ang),0;0,0,1];
H=[1 0 1]; I=eye(3);
for c=1:C
    R=Rtrain(c);
    Q=diag([qOscOverR*R,qOscOverR*R,qBiasOverR*R]);
    dc=fit.dc(c); amp=fit.amplitude(c);
    % Initialise using the training interval; no validation parameter fit.
    x=[amp*sin(fit.phase(c)-2*pi*f0*(trainIdx(1)-1)/fs);...
       -amp*cos(fit.phase(c)-2*pi*f0*(trainIdx(1)-1)/fs);dc];
    % The initial harmonic state is only a rough initialisation; early
    % samples are excluded by the 60%-95% validation window.
    P=diag([max(var(y(trainIdx,c)),1e-9),...
            max(var(y(trainIdx,c)),1e-9),max(R,1e-9)]);
    for n=1:N
        xp=A*x;Pp=A*P*A'+Q;
        innovation=y(n,c)-H*xp;
        kg=(Pp*H')/(H*Pp*H'+R);
        x=xp+kg*innovation;
        P=(I-kg*H)*Pp*(I-kg*H)'+kg*R*kg';
        yhat(n,c)=H*x;
    end
end
end
