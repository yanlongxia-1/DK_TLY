function s=E07D_entropy(modes)
%E07D_ENTROPY 原论文最小包络熵的可重复工程代理，非作者完整适应度源码。
K=size(modes,3);C=size(modes,2);ent=zeros(K,C);
for k=1:K
    for c=1:C
        envelope=abs(gyro_analytic(modes(:,c,k)));
        p=envelope/(sum(envelope)+eps);
        ent(k,c)=-sum(p.*log(p+eps))/log(numel(p));
    end
end
s=mean(min(ent,[],1));
ener=squeeze(sum(sum(modes.^2,1),2));
if any(ener/(sum(ener)+eps)<1e-5),s=s+1;end
if ~isfinite(s),s=Inf;end
end
