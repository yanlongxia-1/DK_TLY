function xout = gyro_haar(x,levels,thresholdScale)
% Orthogonal Haar DWT soft thresholding. Self-contained wavelet baseline.
% Original paper doesn't identify enough wavelet settings for exact cloning.
x=double(x(:)); n=numel(x);
L=min(levels, floor(log2(n)));
if L<1, xout=x; return; end
block=2^L;
N=ceil(n/block)*block;
v=[x; flipud(x(max(1,n-(N-n)+1):n))];
if numel(v)<N
    v=[x; repmat(x(end),N-n,1)];
end
v=v(1:N);
detail=cell(L,1);
for k=1:L
    a=(v(1:2:end)+v(2:2:end))/sqrt(2);
    detail{k}=(v(1:2:end)-v(2:2:end))/sqrt(2);
    v=a;
end
sigma=median(abs(detail{1}-median(detail{1})))/0.67448975;
lambda=thresholdScale*sigma*sqrt(2*log(N));
for k=1:L
    detail{k}=sign(detail{k}).*max(abs(detail{k})-lambda,0);
end
for k=L:-1:1
    vv=zeros(numel(v)*2,1);
    vv(1:2:end)=(v+detail{k})/sqrt(2);
    vv(2:2:end)=(v-detail{k})/sqrt(2);
    v=vv;
end
xout=v(1:n);
end
