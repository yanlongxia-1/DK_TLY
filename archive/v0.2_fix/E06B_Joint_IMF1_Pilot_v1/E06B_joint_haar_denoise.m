function [out, stats] = E06B_joint_haar_denoise(mode, levels, factor)
%E06B_JOINT_HAAR_DENOISE 双通道共享收缩权重的正交 Haar 细节滤波。
% mode: N x 2 的同一联合MVMD模态；只处理第1至levels级高频细节。
% 不使用转台角度/真实运动/测试边带频率；不依赖小波工具箱。
validateattributes(mode,{'numeric'},{'real','finite','2d','ncols',2});
validateattributes(levels,{'numeric'},{'scalar','integer','>=',1,'<=',10});
validateattributes(factor,{'numeric'},{'scalar','nonnegative','finite'});
[N,~] = size(mode);
L = 2^levels;
pad = mod(-N,L);
if pad > 0
    assert(pad < N, '输入过短，无法做镜像填充');
    extended = [mode; flipud(mode(end-pad+1:end,:))];
else
    extended = mode;
end
approx = extended;
details = cell(levels,1);
for j=1:levels
    first = approx(1:2:end,:);
    second = approx(2:2:end,:);
    details{j}=(first-second)/sqrt(2);
    approx=(first+second)/sqrt(2);
end
% 用两路最高频细节的 pooled MAD 估计噪声标准差。
d1=details{1}(:);
sigma = 1.4826*median(abs(d1-median(d1)));
if ~isfinite(sigma),sigma=0;end
numChanged=0;
thresholds=zeros(levels,1);
for j=1:levels
    d = details{j};
    % 同一个时间-尺度位置的a、s共享一个权重，避免两通道单独阈值所致的额外失配。
    threshold = factor*sigma*sqrt(2*log(max(size(d,1)*2,2)));
    thresholds(j)=threshold;
    vectorNorm=sqrt(sum(d.^2,2));
    gain=max(0,1-threshold./max(vectorNorm,eps));
    updated=d.*gain;
    numChanged=numChanged+sum(any(abs(updated-d)>eps,2));
    details{j}=updated;
end
for j=levels:-1:1
    d=details{j};
    restored=zeros(size(approx,1)*2,2);
    restored(1:2:end,:)=(approx+d)/sqrt(2);
    restored(2:2:end,:)=(approx-d)/sqrt(2);
    approx=restored;
end
out=approx(1:N,:);
stats.sigma=sigma;
stats.thresholds=thresholds;
stats.changedPairs=numChanged;
stats.deltaRMS=sqrt(mean((out(:)-mode(:)).^2));
end
