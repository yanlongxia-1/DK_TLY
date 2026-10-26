function yhat=gyro_kf(y,fs,Rtrain,qOscRatio,qBiasRatio)
% Standalone causal oscillator KF comparator, not claimed to match user's older source.
% Each channel contains harmonic in-phase/quadrature states plus slow bias.
N=size(y,1);C=size(y,2); yhat=zeros(N,C);
if ~isscalar(Rtrain) && numel(Rtrain)~=C, error('Bad Rtrain');end
% This implementation uses caller's known/estimated carrier frequency via arg fs(2).
if numel(fs)~=2,error('fs argument must be [sampleRate,carrierHz]');end
sampleRate=fs(1);f0=fs(2);
th=2*pi*f0/sampleRate;
A=[cos(th),-sin(th),0;sin(th),cos(th),0;0,0,1];
H=[1,0,1];I=eye(3);
for c=1:C
    R=max(Rtrain(min(c,numel(Rtrain))),eps);
    Q=diag([qOscRatio*R,qOscRatio*R,qBiasRatio*R]);
    x=[y(1,c)-mean(y(:,c));0;mean(y(:,c))];
    P=diag([max(var(y(:,c)),eps),max(var(y(:,c)),eps),max(var(y(:,c)),eps)]);
    for n=1:N
        xp=A*x;Pp=A*P*A'+Q;
        innovation=y(n,c)-H*xp;
        gain=(Pp*H')/(H*Pp*H'+R);
        x=xp+gain*innovation;
        P=(I-gain*H)*Pp*(I-gain*H)'+gain*R*gain';
        yhat(n,c)=x(1)+x(3);
    end
end
end
