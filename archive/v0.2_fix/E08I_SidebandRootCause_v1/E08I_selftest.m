function E08I_selftest()
cfg=E08I_config(); N=3000; t=(0:N-1)'/cfg.fs;
X=[.1+cos(2*pi*8*t)+.3*cos(2*pi*9*t)+.1*cos(2*pi*10*t), ...
   -.1+sin(2*pi*8*t)+.3*sin(2*pi*9*t)+.1*sin(2*pi*10*t)];
m=E07B_candidate_metrics(X,X,t,[8 9 10],[8 9 10],cfg);
assert(max(abs(m.retentionPct(:)-100))<1e-7 && m.maxPhaseShiftDeg<1e-7 && m.maxRelativePhaseShiftDeg<1e-7);
assert(m.sidebandPass, '恒等重构物理保真自检失败');
% 不声称验证MVMD实际收敛，仅验证调用接口和形状
[U,c,info,~]=E08G_mvmd_damped(X,cfg.fs,3,154.157959862,25,cfg.tol,cfg.tau,9,1,25);
assert(isequal(size(U),[N 2 3]) && numel(c)==3 && isfield(info,'ACreconRelative'));
fprintf('PASS E08I: 物理指标与求解器函数调用自检通过；不代表真实实验已通过。\n');
end
