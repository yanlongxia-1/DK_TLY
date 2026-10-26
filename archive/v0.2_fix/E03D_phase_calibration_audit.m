%% E03D - Dual-orthogonal gyro: phase / bias / sensitivity error attribution
% Prerequisites: v0.3 MATLAB project, gyro_config.m, gyro_mvmd.m,
% gyro_reconstruct.m, and converged E02 results/record.mat.
%
% Scientific guardrails:
%  (1) K, alpha and IMF thresholds FROZEN from E02 (no tuning/optimization).
%  (2) New fixed simulation seeds not used by E03B/E03C.
%  (3) 'oracle' knows true phi, bias and gains ONLY as a theoretical reference.
%      The 'nominal' decoder assumes fixed 3 Hz phase, original bias and gains.
%  (4) Filter may suppress noise, but cannot infer unknown phi from two channels
%      without additional assumptions/measurements.
%  (5) NOT a full physical-constraint optimization or a real-sensor accuracy test.
%  (6) At most 8 MVMD runs, no FVIM, immediate checkpoint after each case.
clear; clc; close all;
cfgNow=gyro_config();
[~,stem]=fileparts(cfgNow.file);
recordPath=fullfile(cfgNow.outputRoot,stem,'E02','record.mat');
if ~isfile(recordPath)
    error('E03D:noRecord','No E02 record.mat: %s',recordPath);
end
S=load(recordPath,'record'); r=S.record;
if ~isfield(r,'mvmdInfo') || ~r.mvmdInfo.converged
    error('E03D:baseline','A converged E02 baseline is required.');
end
cfg=r.cfg; fs=cfg.fs; K=r.k; alpha=r.alpha;
if abs(fs-1000)>1e-9
    error('E03D:fs','This preregistered experiment assumes 1000 Hz.');
end
N=2048; t=(0:N-1)'/fs; T=t(end);
Ka0=0.0020; Ks0=0.0022;
b0=[-0.186041, -0.072559];
froll=3.0; phiNom=2*pi*froll*t;
ox=20+8*sin(2*pi*3*t);
oy=8+6*cos(2*pi*3*t);
sigma=0.012;
scenarios={'ideal','phase_offset_5deg','roll_drift_0p10Hz','bias_gain_mismatch'};
seeds=[20271001,20271002];
outputDir=fullfile(cfgNow.outputRoot,stem,'E03D');
if ~isfolder(outputDir), mkdir(outputDir); end
resultFile=fullfile(outputDir,'E03D_phase_calibration_results.csv');
notesFile=fullfile(outputDir,'E03D_notes.txt');
if isfile(resultFile)
    error('E03D:exists',['Output already exists. Rename the E03D directory ', ...
        'to avoid mixing runs before starting again.']);
end
cols={'scenario','seed','method','converged','iterations','last_delta', ...
    'AC_recon_relative','runtime_s','phi_end_error_deg','bias_error_a_units', ...
    'bias_error_s_units','gain_error_a_pct','gain_error_s_pct', ...
    'OmegaX_RMSE_deg_s','OmegaY_RMSE_deg_s','OmegaMagnitude_RMSE_deg_s', ...
    'nClassChanges','nPhysicsProtected','max_E02_E03_sensor_delta'};
rows=cell(0,numel(cols));
clockAll=tic; caseNo=0;
for ci=1:numel(scenarios)
    scenario=scenarios{ci};
    for si=1:numel(seeds)
        seed=seeds(si);
        rng(seed,'twister');
        phiTrue=phiNom; Ktrue=[Ka0,Ks0]; bTrue=b0;
        switch scenario
            case 'ideal'
                % All assumed calibration parameters happen to be correct.
            case 'phase_offset_5deg'
                phiTrue=phiNom+5*pi/180;
            case 'roll_drift_0p10Hz'
                % f_roll(t)=3+0.10*t/T; final frequency rises to 3.1 Hz.
                phiTrue=2*pi*(froll*t+0.5*0.10*t.^2/T);
            case 'bias_gain_mismatch'
                Ktrue=[Ka0*1.06,Ks0*0.94];
                bTrue=b0+[0.004,-0.003];
        end
        clean=[Ktrue(1)*(ox.*sin(phiTrue)-oy.*cos(phiTrue)), ...
               Ktrue(2)*(ox.*cos(phiTrue)+oy.*sin(phiTrue))];
        measured=clean+bTrue+sigma*randn(N,2);
        % Identical impulse locations/amplitudes for every paired case.
        ix=[402,1101,1719];
        measured(ix,1)=measured(ix,1)+[0.055;-0.047;0.050];
        measured(ix+9,2)=measured(ix+9,2)+[-0.052;0.041;-0.056];
        phiEndError=(phiTrue(end)-phiNom(end))*180/pi;
        tStart=tic;
        [modes,centers,info]=gyro_mvmd(measured,fs,K,alpha, ...
            min(1100,cfg.mvmd.finalIterations),cfg.mvmd.tolerance,cfg.mvmd.tau);
        elapsed=toc(tStart);
        filt02=nan(N,2); filt03=nan(N,2);
        changed=NaN; protected=NaN; delta03=NaN;
        if info.converged
            [filt02,detail02]=gyro_reconstruct(modes,centers,measured,cfg,'paper');
            [filt03,detail03]=gyro_reconstruct(modes,centers,measured,cfg,'physics');
            changed=sum(detail02.class_0drop_1mix_2keep ~= ...
                        detail03.class_0drop_1mix_2keep);
            protected=sum(detail03.physical_protection_applied);
            delta03=max(abs(filt03(:)-filt02(:)));
        else
            warning('E03D:notConverged','%s seed=%d MVMD stopped at %d iterations', ...
                    scenario,seed,info.iterations);
        end
        dataMethods={measured,measured,filt02,filt02,filt02,filt03};
        labels={'raw_oracle','raw_nominal','E02_oracle', ...
            'E02_nominal_phi','E02_nominal_all','E03_nominal_all'};
        usedPhi={phiTrue,phiNom,phiTrue,phiNom,phiNom,phiNom};
        usedB={bTrue,b0,bTrue,bTrue,b0,b0};
        usedK={Ktrue,[Ka0,Ks0],Ktrue,Ktrue,[Ka0,Ks0],[Ka0,Ks0]};
        roi=(101:N-100)';
        for mi=1:numel(labels)
            Y=dataMethods{mi};
            bx=usedB{mi}; kx=usedK{mi}; ph=usedPhi{mi};
            if all(isfinite(Y(:)))
                a=(Y(:,1)-bx(1))/kx(1); s=(Y(:,2)-bx(2))/kx(2);
                xhat=s.*cos(ph)+a.*sin(ph);
                yhat=s.*sin(ph)-a.*cos(ph);
                errX=sqrt(mean((xhat(roi)-ox(roi)).^2));
                errY=sqrt(mean((yhat(roi)-oy(roi)).^2));
                errMag=sqrt(mean((hypot(xhat(roi),yhat(roi))-hypot(ox(roi),oy(roi))).^2));
            else
                errX=NaN; errY=NaN; errMag=NaN;
            end
            rows(end+1,:)={scenario,seed,labels{mi},double(info.converged), ...
                info.iterations,info.delta,info.AC_reconRelative,elapsed, ...
                phiEndError,bTrue(1)-b0(1),bTrue(2)-b0(2), ...
                100*(Ktrue(1)/Ka0-1),100*(Ktrue(2)/Ks0-1), ...
                errX,errY,errMag,changed,protected,delta03}; %#ok<AGROW>
        end
        Q=cell2table(rows,'VariableNames',cols);
        writetable(Q,resultFile);  % Checkpoint after every case
        caseNo=caseNo+1;
        fprintf('E03D %d/8 %-23s seed=%d converged=%d MVMD=%.2fs\n', ...
            caseNo,scenario,seed,info.converged,elapsed);
        clear modes centers detail02 detail03
        if toc(clockAll)>120 && caseNo<8
            warning('E03D:softLimit','Soft runtime cap 120s reached. Saved completed cases.');
            break;
        end
    end
    if toc(clockAll)>120, break; end
end
fid=fopen(notesFile,'w');
fprintf(fid,'E03D: phase/bias/gain uncertainty sensitivity; 4 scenarios x 2 held-out seeds.\n');
fprintf(fid,'Fixed trained E02 K=%d alpha=%.10f, tau=%g, maximum 1100 iterations.\n',K,alpha,cfg.mvmd.tau);
fprintf(fid,'No FVIM, no tuning, no use of E03C simulation truth.\n');
fprintf(fid,'Oracle knows true phase, sensitivities, and biases ONLY for error attribution.\n');
fprintf(fid,'Nominal assumes 3 Hz and zero initial phase and nominal calibration.\n');
fprintf(fid,'Scenario roll drift is 3 to 3.10 Hz over 2.047 seconds.\n');
fprintf(fid,'One or more nonconverged cases must remain in convergence statistics.\n');
fprintf(fid,'RMSE units are deg/s, synthetic truth only, ROI samples 101:1948.\n');
fprintf(fid,'Interpretation: E02_oracle vs raw_oracle measures denoising potential.\n');
fprintf(fid,'E02_nominal_phi vs E02_oracle isolates roll phase error.\n');
fprintf(fid,'E02_nominal_all vs E02_nominal_phi adds calibration mismatch.\n');
fprintf(fid,'E03_nominal_all equals E02_nominal_all if guard never activates.\n');
fprintf(fid,'Code static-reviewed, not MATLAB-runtime tested.\n');
fclose(fid);
fprintf('E03D complete: %s\n',resultFile);
