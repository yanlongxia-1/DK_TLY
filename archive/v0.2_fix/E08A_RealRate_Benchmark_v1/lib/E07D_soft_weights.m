function [keepWeight,mixWeight]=E07D_soft_weights(cs,samp,cfg)
%E07D_SOFT_WEIGHTS CS/SampEn双通道共享、连续可导的模式保留/混叠权重。
% 共享权重仅根据两路联合的cs和samp，不使用转台真值。
% 固定门槛中心与原论文工程版相同；过渡宽度是新增超参数。
cs=double(cs);samp=double(samp);
assert(isequal(size(cs),size(samp)),'cs与samp的大小应一致');
assert(all(isfinite(cs(:))) && all(isfinite(samp(:))),'模态指标不允许NaN');
sgm=@(z) 1./(1+exp(-max(-50,min(50,z))));
keepWeight=sgm((cs-cfg.class.csKeep)/cfg.soft.csWidth) .* ...
           sgm((cfg.class.sampKeep-samp)/cfg.soft.sampWidth);
mixGate=sgm((cs-cfg.class.csMixed)/cfg.soft.csWidth) .* ...
        sgm((cfg.class.sampMixed-samp)/cfg.soft.sampWidth);
mixWeight=(1-keepWeight).*mixGate.^cfg.soft.mixExponent;
assert(all(keepWeight(:)>=0 & keepWeight(:)<=1) && ...
       all(mixWeight(:)>=0 & mixWeight(:)<=1),'权重超出物理范围');
end
