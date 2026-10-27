function E08H_selftest()
% Only source-level data flow and structure test; does not imply actual convergence.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
fs=1000; t=(0:1999)'/fs;
X=[.04*cos(2*pi*9*t)+.001*sin(2*pi*57*t),...
   .04*sin(2*pi*9*t)+.001*cos(2*pi*47*t)];
[U,center,info,snaps]=E08G_mvmd_damped(X,fs,3,154,35,1e-5,0,9,1,35);
assert(isequal(size(U),[2000 2 3]) && numel(center)==3);
assert(~isempty(snaps) && all(isfinite(U(:))) && isfinite(info.ACreconRelative));
fprintf('E08H selftest PASS: MATLAB solver and snapshot fields accessible.\n');
end
