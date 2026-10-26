function out = gyro_sinefit(x, fs, f0)
% Harmonic least-squares including DC; phase convention: A*sin(2*pi*f*t+phi)
N=size(x,1); t=(0:N-1)'/fs;
H=[sin(2*pi*f0*t),cos(2*pi*f0*t),ones(N,1)];
B=H\double(x);
out.amplitude=hypot(B(1,:),B(2,:));
out.phase=atan2(B(2,:),B(1,:));
out.dc=B(3,:);
out.predicted=H*B;
out.residual=double(x)-out.predicted;
out.snrProxy=10*log10(max(sum((out.predicted-mean(out.predicted,1)).^2,1),eps)./max(sum(out.residual.^2,1),eps));
end
