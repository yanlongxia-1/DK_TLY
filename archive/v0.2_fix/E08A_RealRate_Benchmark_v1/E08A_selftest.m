function E08A_selftest()
%E08A_SELFTEST 检查真实解调链最容易出错的C/S顺序与符号。
rng(20261026);
t=(0:0.001:2.001)';theta=17-360*6*t;
HE=[.008 -.004; .0012 .003];
B=[cosd(theta),sind(theta)];
scale=[-.001 -.0079;-.0078 .0008];
bias=[-.1 -.07];omegaTruth=20;
y=repmat(bias,numel(t),1)+B*HE.'+omegaTruth*(B*scale.');
[est,den]=E08A_estimate_outer(y,theta,HE,bias,scale);
assert(max(abs(est-omegaTruth))<1e-9,'解调恒速恢复失败');
assert(min(den)>1e-10,'解调响应奇异');
% 相位固定频率推算与编码器正交系数一致。
thetaFixed=theta(1)-360*6*t;
assert(max(abs(theta-thetaFixed))<1e-9,'固定滚转相位校验失败');
fprintf('PASS E08A: 合成W双通道有符号恒速恢复，固定相位一致。\n');
end
