function [omega,den]=E08A_estimate_outer(Y,thetaDeg,HE,bias,scale)
%E08A_ESTIMATE_OUTER E类零输入基线+双路复响应的W恒速等效输入估计。
% 注意 scale = [Ca Sa; Cs Ss]，H=C+iS仅为系数记法，不能用复数直接相除。
assert(size(Y,2)==2 && size(HE,1)==2 && size(HE,2)==2,...
    '输入必须为双通道及2x2谐波系数');
assert(numel(bias)==2 && numel(scale)==4,'标定维数错误');
thetaDeg=thetaDeg(:);
assert(size(Y,1)==numel(thetaDeg),'相位与信号长度不一致');
B=[cosd(thetaDeg),sind(thetaDeg)];
g=B*scale.';
u=Y-repmat(reshape(bias,1,2),numel(thetaDeg),1)-B*HE.';
den=sum(g.^2,2);
assert(all(isfinite(den)) && min(den)>1e-12,'响应分母接近0，无法安全解调');
omega=sum(g.*u,2)./den;
end
