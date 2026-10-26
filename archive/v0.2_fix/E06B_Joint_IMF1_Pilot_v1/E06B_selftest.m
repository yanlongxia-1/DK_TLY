function E06B_selftest()
%E06B_SELFTEST 检查Haar分解+逆变换无阈值时精确恢复双路信号。
n=4000; t=(0:n-1)'/1000;
x=[sin(2*pi*2*t)+0.3*cos(2*pi*4*t),cos(2*pi*2*t)+0.2*sin(2*pi*4*t)];
[y,~]=E06B_joint_haar_denoise(x,6,0);
assert(max(abs(x(:)-y(:)))<2e-12, 'E06B Haar零阈值可逆性自检失败');
fprintf('PASS: 双通道Haar零阈值自检。\n');
end
