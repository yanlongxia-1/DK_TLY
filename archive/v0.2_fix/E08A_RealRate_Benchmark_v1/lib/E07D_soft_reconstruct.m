function [Y,details]=E07D_soft_reconstruct(modes,centerHz,X,cfg)
%E07D_SOFT_RECONSTRUCT 固定MVMD模态，仅连续改变联合IMF重构。
% 共享权重作用于双通道的成对IMF；混叠小波基础处理保留原工程行为。
% 无转台角速度/姿态真值进入本函数。
[N,C,K]=size(modes);
assert(C==2 && isequal(size(X),[N C]),'输入必须为N×2');
[~,details]=gyro_reconstruct(modes,centerHz,X,cfg,'paper');
[keep,mix]=E07D_soft_weights(details.cosine_abs,details.sampEn,cfg);
Y=repmat(mean(X,1),N,1);
for k=1:K
    M=modes(:,:,k);
    W=zeros(N,C);
    for c=1:C
        W(:,c)=gyro_haar(M(:,c),cfg.wavelet.levels,cfg.wavelet.thresholdScale);
    end
    Y=Y+keep(k)*M+mix(k)*W;
end
details.softKeepWeight=keep;
details.softMixWeight=mix;
assert(all(isfinite(Y(:))),'软重构输出出现非有限值');
end
