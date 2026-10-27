function E08G_selftest()
%E08G_SELFTEST Quick synthetic path check, no user TXT files or lengthy FVIM search.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'lib'));
cfg=E08G_config();
fs=cfg.fs;
t=(0:999)'/fs;
X=[0.15*cos(2*pi*cfg.rollHz*t)+0.003*cos(2*pi*27*t), ...
   0.15*sin(2*pi*cfg.rollHz*t)+0.003*sin(2*pi*27*t)];
for eta=[1.0, 0.5]
  [modes,centerHz,info,snaps]=E08G_mvmd_damped(X,fs,3,154.157959862,...
     35,1e-5,0,cfg.rollHz,eta,[10 35]);
  assert(~isempty(snaps) && isfield(snaps,'centerStepHz'), ...
     'Snapshot centerStepHz is missing');
  assert(isfield(info,'ACreconRelative') && isfinite(info.ACreconRelative), ...
     'Solver output ACreconRelative is missing/nonfinite');
  assert(isfield(info,'centerStepHz') && isfinite(info.centerStepHz), ...
     'Solver output centerStepHz is missing/nonfinite');
  assert(size(modes,1)==size(X,1) && size(modes,2)==2, ...
     'Unexpected decomposed modes shape');
  [Y,details]=gyro_reconstruct(snaps(end).modes,centerHz,X,cfg,'paper');
  assert(all(isfinite(Y(:))) && height(details)==3, ...
     'Reconstruction/checkpoint path failed');
  fprintf('PASS eta=%.1f; iter=%d; converged=%d; ACerr=%.4g; snaps=%d\n',...
     eta,info.iterations,info.converged,info.ACreconRelative,numel(snaps));
end
fprintf('PASS E08G_selftest: struct fields / snapshot reconstruction are consistent.\n');
end
