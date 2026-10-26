%% E03C - Dynamic sideband / IMF occupancy audit (diagnostic, low computational load)
% Research objective:
%   Locate real transverse-rate modulation information in paired MVMD modes.
%   Do NOT tune K, alpha, or thresholds against this diagnostic set.
%   This does not estimate roll phase or bias from real data.
% Run from the v0.3 MATLAB project root with gyro_config.m and the
% gyro_mvmd.m / gyro_reconstruct.m functions on MATLAB path.
%
% Safety limits: 8 decompositions total, no FVIM, 2048 samples/case,
% max 800 iterations, write outputs after every completed case.
clear; clc; close all;
cfgNow=gyro_config();
[~, stem]=fileparts(cfgNow.file);
recordFile=fullfile(cfgNow.outputRoot,stem,'E02','record.mat');
if ~isfile(recordFile)
    error('E03C:noRecord','E02 record.mat missing: %s',recordFile);
end
S=load(recordFile,'record'); r=S.record;
if ~isfield(r,'mvmdInfo') || ~r.mvmdInfo.converged
    error('E03C:invalidBaseline','Saved E02 must have converged.');
end
cfg=r.cfg; fs=cfg.fs;
if abs(fs-1000)>1e-9
    error('E03C:fsMismatch','This preset assumes fs=1000 Hz; inspect before changing.');
end
K=r.k; alpha=r.alpha;
N=2048; t=(0:N-1)'/fs;
Ka=0.0020; Ks=0.0022;
bias=[-0.186041,-0.072559]; % Synthetic-only known bias
sigma=0.012;
modFreqs=[1,3,8,15]; % registered diagnostic probes, NOT an optimized frequency search
seeds=[20261111,20261112];
outputDir=fullfile(cfgNow.outputRoot,stem,'E03C');
if ~isfolder(outputDir),mkdir(outputDir);end
perfFile=fullfile(outputDir,'E03C_true_rate_diagnostics.csv');
modeFile=fullfile(outputDir,'E03C_mode_occupancy.csv');
if isfile(perfFile) || isfile(modeFile)
    error('E03C:exists',['E03C output exists; to prevent mixed runs, rename the E03C folder ',...
        'before starting a new run.']);
end
perfRows={}; modeRows={};
completed=0;
for fi=1:numel(modFreqs)
  fmod=modFreqs(fi);
  for si=1:numel(seeds)
    seed=seeds(si);
    rng(seed+100*fi,'twister');
    phi=2*pi*3*t; % true phi is known ONLY in this synthetic experiment
    ox=20+8*sin(2*pi*fmod*t);
    oy=8+6*cos(2*pi*fmod*t);
    clean=[Ka*(ox.*sin(phi)-oy.*cos(phi)), ...
           Ks*(ox.*cos(phi)+oy.*sin(phi))];
    cleanDC=[Ka*(20*sin(phi)-8*cos(phi)), ...
             Ks*(20*cos(phi)+8*sin(phi))];
    sidebandClean=clean-cleanDC;
    meas=clean+bias+sigma*randn(N,2);
    spikes=[402,1101,1719];
    meas(spikes,1)=meas(spikes,1)+[0.055;-0.047;0.05];
    meas(spikes+9,2)=meas(spikes+9,2)+[-0.052;0.041;-0.056];
    tic;
    [modes, centers, info]=gyro_mvmd(meas,fs,K,alpha, ...
        min(800,cfg.mvmd.finalIterations),cfg.mvmd.tolerance,cfg.mvmd.tau);
    elapsed=toc;
    if info.converged
       [filter02, detail02]=gyro_reconstruct(modes,centers,meas,cfg,'paper');
       [filter03, detail03]=gyro_reconstruct(modes,centers,meas,cfg,'physics');
    else
       filter02=nan(size(meas)); filter03=nan(size(meas));
       detail02=[]; detail03=[];
       warning('E03C:notConverged','fmod=%g seed=%d not converged: diagnostics invalid',fmod,seed);
    end
    signals={meas,filter02,filter03}; labels={'raw','E02_fixed','E03_fixed_physics'};
    ii=(101:N-100)'; % identical boundary trimming for all methods
    for methodID=1:numel(signals)
       ys=signals{methodID};
       if all(isfinite(ys(:)))
           aa=(ys(:,1)-bias(1))/Ka; ss=(ys(:,2)-bias(2))/Ks;
           xhat=ss.*cos(phi)+aa.*sin(phi);
           yhat=ss.*sin(phi)-aa.*cos(phi);
           eX=sqrt(mean((xhat(ii)-ox(ii)).^2));
           eY=sqrt(mean((yhat(ii)-oy(ii)).^2));
           eMag=sqrt(mean((hypot(xhat(ii),yhat(ii))-hypot(ox(ii),oy(ii))).^2));
           % Dynamic gain: least-squares coefficient at known modulation frequency.
           H=[ones(numel(ii),1),sin(2*pi*fmod*t(ii)),cos(2*pi*fmod*t(ii))];
           cx=H\xhat(ii);cy=H\yhat(ii);
           gainX=hypot(cx(2),cx(3))/8;
           gainY=hypot(cy(2),cy(3))/6;
           [~,ixBest]=min(abs(centers-3));
           if info.converged
               changed=sum(detail02.class_0drop_1mix_2keep~=detail03.class_0drop_1mix_2keep);
               protected=sum(detail03.physical_protection_applied);
               outputDelta=max(abs(filter03(:)-filter02(:)));
           else
               changed=NaN;protected=NaN;outputDelta=NaN;
           end
       else
           eX=NaN;eY=NaN;eMag=NaN;gainX=NaN;gainY=NaN;
           changed=NaN;protected=NaN;outputDelta=NaN;ixBest=NaN;
       end
       perfRows(end+1,:)={fmod,seed,labels{methodID},K,alpha,double(info.converged),...
           info.iterations,elapsed,info.AC_reconRelative,eX,eY,eMag,gainX,gainY,...
           changed,protected,outputDelta}; %#ok<AGROW>
    end
    if info.converged
       sidePow=sum(sidebandClean(:).^2)+eps;
       for k=1:K
          mk=modes(:,:,k);
          % Signed normalized projection is DIAGNOSTIC using simulation truth.
          % Not a usable real-world classification statistic.
          projectedSideband=sum(mk(:).*sidebandClean(:))/sidePow;
          modeRows(end+1,:)={fmod,seed,k,centers(k),...
             detail02.class_0drop_1mix_2keep(k),...
             detail03.class_0drop_1mix_2keep(k),...
             detail03.physical_protection_applied(k),...
             detail02.energy_share(k), projectedSideband,...
             detail03.physics_band_share(k)}; %#ok<AGROW>
       end
    end
    completed=completed+1;
    P=cell2table(perfRows,'VariableNames',{'modulation_Hz','seed','method','K_fixed',...
      'alpha_fixed','converged','iterations','runtime_s','AC_recon_relative',...
      'TRUE_OmegaX_RMSE','TRUE_OmegaY_RMSE','TRUE_OmegaMagnitude_RMSE',...
      'OmegaX_modulation_gain','OmegaY_modulation_gain',...
      'nChangedClasses','nProtectedModes','max_E02_E03_delta'});
    writetable(P,perfFile); % checkpoint immediately
    if ~isempty(modeRows)
       M=cell2table(modeRows,'VariableNames',{'modulation_Hz','seed','IMF',...
         'center_Hz','E02_class','E03_class','physical_protection',...
         'mode_energy_share','true_sideband_projection_diagnostic','physics_band_share'});
       writetable(M,modeFile);
    end
    fprintf('E03C: %d/8; modulation=%.1fHz; seed=%d; converged=%d; %.2fs\n',...
       completed,fmod,seed,info.converged,elapsed);
  end
end
fid=fopen(fullfile(outputDir,'E03C_notes.txt'),'w');
fprintf(fid,'E03C: Fixed real-data-trained K, alpha, classifier thresholds; no FVIM.\n');
fprintf(fid,'Synthetic diagnostics only: known phi, Ka, Ks, biases, and true signed rates.\n');
fprintf(fid,'Modulation frequencies [1,3,8,15]Hz, seeds [20261111,20261112].\n');
fprintf(fid,'Measured RMSE is computed on indices 101:%d to reduce boundary effects.\n',N-100);
fprintf(fid,'TRUE sideband projection uses unavailable-in-practice clean simulated signal. DO NOT optimize from it.\n');
fprintf(fid,'This test is EXPLORATORY; create fresh randomized holdout before final claims.\n');
fprintf(fid,'Mode classes: 0=discard, 1=wavelet-denoise, 2=keep.\n');
fprintf(fid,'One-to-one comparison E02/E03 is meaningful ONLY if the physics guard actually changes output.\n');
fprintf(fid,'Do not interpret numeric solver convergence as accuracy or real-time readiness.\n');
fclose(fid);
fprintf('E03C DONE. Return E03C_true_rate_diagnostics.csv, E03C_mode_occupancy.csv, E03C_notes.txt\n');
