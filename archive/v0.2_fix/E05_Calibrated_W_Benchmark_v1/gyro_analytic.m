function z = gyro_analytic(x)
% Hilbert analytic transform implemented by FFT; no toolbox needed.
N = size(x,1);
X = fft(x, [], 1);
h = zeros(N,1);
if mod(N,2)==0
    h([1,N/2+1])=1;
    h(2:N/2)=2;
else
    h(1)=1;
    h(2:(N+1)/2)=2;
end
z = ifft(X.*h, [], 1);
end
