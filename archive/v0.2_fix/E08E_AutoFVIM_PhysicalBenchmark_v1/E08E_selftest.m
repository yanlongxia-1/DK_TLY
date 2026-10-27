function E08E_selftest()
%E08E_SELFTEST 合成物理解调及数值算法小样检查（非真实精度验证）。
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
fs=1000;t=(0:1199)'/fs;f=9;
X=[.02*cos(2*pi*f*t)+.001*sin(2*pi*90*t), .02*sin(2*pi*f*t)+.001*cos(2*pi*87*t)];
[m,omega,info]=gyro_mvmd_conditioned(X,fs,3,180,40,1e-5,0,f,true);
assert(isequal(size(m),[numel(t),2,3])&&all(isfinite(m(:)))&&all(isfinite(omega(:))));
assert(info.conditionedInitialization&&any(abs(info.initialCenterHz-f)<1e-12));
theta=17-360*f*t;B=[cosd(theta) sind(theta)];
scale=[-.001 -.0079;-.0078 .0008];HE=[.008 -.004; .0012 .003];b=[-.1 -.07];
Y=repmat(b,numel(t),1)+B*HE.'+20*(B*scale.');
[rate,den]=E08A_estimate_outer(Y,theta,HE,b,scale);
assert(max(abs(rate-20))<1e-9&&min(den)>1e-10);
assert(isfinite(E07D_entropy(m)));
fprintf('PASS E08E: 数值维度、频率初始化、熵目标及有符号解调自检通过。\n');
end
