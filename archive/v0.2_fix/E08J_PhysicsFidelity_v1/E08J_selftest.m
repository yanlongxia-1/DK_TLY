function E08J_selftest()
% Fast, no real data or reference encoder needed.
cfg=E08J_config();N=3000;t=(0:N-1)'/cfg.fs;f=[8 9 10];
X=[.040*cos(2*pi*8*t)+.001*cos(2*pi*9*t)+.037*cos(2*pi*10*t),...
   .039*sin(2*pi*8*t)+.0012*sin(2*pi*9*t)+.038*sin(2*pi*10*t)];
% Original 1% ratio fails for weak 9 Hz while new ABSOLUTE constraint passes.
Y=X + [.00004*cos(2*pi*9*t),.000045*sin(2*pi*9*t)];
[o,detail]=E08J_physics_metrics(X,Y,t,f,cfg);
assert(~o.oldPass && o.newPass,'弱分量判据回归失败');
assert(height(detail)==6 && all(isfinite(detail.complexCoefficientError)));
% 2% damage to a physically strong sideband MUST still be rejected.
Y2=X + [.0008*cos(2*pi*8*t),zeros(N,1)];
[o2,~]=E08J_physics_metrics(X,Y2,t,f,cfg);
assert(~o2.newPass,'强边带损伤不应被新规则放行');
% Strong-component phase loss must also be rejected.
Y3=[.040*cos(2*pi*8*t+.02)+.001*cos(2*pi*9*t)+.037*cos(2*pi*10*t),X(:,2)];
[o3,~]=E08J_physics_metrics(X,Y3,t,f,cfg);
assert(~o3.newPass,'强边带相位漂移不应被放行');
% Check integration with existing strictly configured MVMD engine.
Xshort=X(1:300,:);[U,cent,inf,~]=E08G_mvmd_damped(Xshort,cfg.fs,3,154.157959862,40,1e-5,0,9,1,40);
assert(isequal(size(U),[300 2 3]) && all(isfinite(U(:))) && isfield(inf,'ACreconRelative'));
[y,imf]=gyro_reconstruct(U,cent,Xshort,cfg,'paper');
assert(all(isfinite(y(:))) && height(imf)==3);
disp('PASS E08-J selftest: weak absolute fidelity / strong amplitude / strong phase / MVMD chain');
end
