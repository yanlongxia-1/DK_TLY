function E09A0_selftest()
%E09A0_SELFTEST 不涉及外部TXT文件，验证时间索引、单位与积分常值误差。
cfg=E09A0_config();cfg.fs=1000;cfg.sampleTimeToleranceSec=1e-7;
N=5000;t=(0:N-1)'/cfg.fs;
refAngle=13 + 60*t; % 真实外框60deg/s，初始13deg
ix=1001:4000;ts=t(ix);M=numel(ix);
rates=table(ts,60*ones(M,1),60*ones(M,1),61*ones(M,1), ...
  'VariableNames',{'time_s','reference_outer_deg_s','correct','plus1bias'});
[T,X]=E09A0_calculate(refAngle,rates,cfg);
a=T(T.method=="correct",:);b=T(T.method=="plus1bias",:);
assert(height(T)==2 && a.RMSE_encoder_deg < 1e-9);
assert(abs(b.EndError_encoder_deg-(M-1)/cfg.fs)<1e-8);
assert(all(X.originalSampleIndex==ix(:)));
fprintf('PASS E09A0: TXT索引/deg-s积分/编码器角度RMSE/固定偏置末端累计误差。\n');
end
