function E08K_selftest()
% Fast functional checks; not proof of real-world algorithm accuracy.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
cfg = E08K_config();

% 自检使用9Hz滚转频率
cfg.rollHz = cfg.frequenciesHz(1);

N = 3000;
t = (0:N-1)' / cfg.fs;
f=[8 9 10];
X=[.040*cos(2*pi*8*t)+.001*cos(2*pi*9*t)+.037*cos(2*pi*10*t),...
   .039*sin(2*pi*8*t)+.0012*sin(2*pi*9*t)+.038*sin(2*pi*10*t)];
Y=X+[.00004*cos(2*pi*9*t),.000045*sin(2*pi*9*t)];
[q,~]=E08J_physics_metrics(X,Y,t,f,cfg);
assert(~q.oldPass && q.newPass,'Weak-frequency acceptance test failed');
Y2=X+[.0008*cos(2*pi*8*t),zeros(N,1)];
[q2,~]=E08J_physics_metrics(X,Y2,t,f,cfg);
assert(~q2.newPass,'Damaged strong sideband must fail');
[U,cent,info]=E08G_mvmd_damped(X(1:300,:),cfg.fs,3,180,...
 40,1e-5,0,9,1,[]);
assert(isequal(size(U),[300 2 3])&&all(isfinite(U(:)))&&isfield(info,'ACreconRelative'));
[rec,~]=gyro_reconstruct(U,cent,X(1:300,:),cfg,'paper');
assert(all(isfinite(rec(:)))&&isfinite(E07D_entropy(U)));
theta=17-360*9*(0:1199)'/cfg.fs;B=[cosd(theta) sind(theta)];
scale=[-.001 -.0079;-.0078 .0008];HE=[.008 -.004;.0012 .003];bias=[-.1 -.07];
Z=repmat(bias,numel(theta),1)+B*HE.'+20*(B*scale.');
[rate,den]=E08A_estimate_outer(Z,theta,HE,bias,scale);
assert(max(abs(rate-20))<1e-9&&min(den)>1e-10);
fprintf('PASS E08-K: fixed fidelity policy / MVMD / hard reconstruction / entropy / W-equivalent-rate chain.\n');
end
