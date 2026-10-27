function E08C_selftest()
% 数值自检：频率条件化初始化没有改变MVMD输出维度/有限性。
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
fs=1000;t=(0:1199)'/fs;f=9;
X=[.02*cos(2*pi*f*t)+.001*sin(2*pi*90*t),... 
   .02*sin(2*pi*f*t)+.001*cos(2*pi*87*t)];
[m,w,info]=gyro_mvmd_conditioned(X,fs,3,180,40,1e-5,0,9,true);
assert(isequal(size(m),[numel(t),2,3]) && all(isfinite(m(:))) && all(isfinite(w(:))));
assert(info.conditionedInitialization && any(abs(info.initialCenterHz-9)<1e-12));
fprintf('PASS E08C: 已知频率引导初始化、双通道输出维度及数值有限性。\n');
end
