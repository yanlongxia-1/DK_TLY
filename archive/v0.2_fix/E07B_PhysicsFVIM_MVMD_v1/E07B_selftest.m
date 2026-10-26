function E07B_selftest()
% 轻量合成单元检查：不运行MVMD，不读取Z3-2，不声称真实滤波结果。
cfg=E07B_config();fs=cfg.fs;
t=(0:1/fs:3-1/fs)';
X=[0.15*cos(4*pi*t)+0.21*cos(8*pi*t), 0.16*sin(4*pi*t)+0.2*sin(8*pi*t)];
Y=X;
m=E07B_candidate_metrics(X,Y,t,[2 3 4],[2 4],cfg);
assert(all(abs(m.retentionPct(:)-100)<1e-6),'自检:原样输入不应改变幅值');
assert(m.maxPhaseShiftDeg<1e-6,'自检:原样输入不应改变相位');
assert(all(abs(m.nonTargetSuppressionDb)<1e-6),'自检:原样输入不应改变残差代理');
Y2=1.05*X;
m2=E07B_candidate_metrics(X,Y2,t,[2 3 4],[2 4],cfg);
assert(~m2.sidebandPass && m2.maxAmplitudeDeviationPct>4.9,'自检:幅值失真应触发保护');
assert(cfg.search.maxUniqueTrainSolves==12,'自检:12次训练预算已改变');
fprintf('E07-B synthetic score/selftest: PASS\n');
end
