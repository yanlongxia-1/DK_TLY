%% E04: REAL-DATA independent roll-reference feasibility audit (no filtering)
% Put this file beside gyro_config.m and gyro_read.m, and set cfg.file to the
% genuine 7-column 20200824 TXT. Requires MATLAB R2018b or newer.
% NO FVIM/MVMD, NO calibration assumption, NO synthetic truth.
% Col 1-3: outer/middle/inner gimbal angular rates (reported deg/s).
% Col 4: a gyro, Col 6: one accelerometer channel, Col 7: s gyro.
% The audit does NOT claim col3 equals body roll rate or col6 is roll angle.
clear; clc; close all;
cfg = gyro_config();
d = gyro_read(cfg.file,cfg);
A = double(d.A);
if size(A,2)<7, error('E04 requires at least 7 numeric columns.'); end
if any(~isfinite(A(:,[1 2 3 4 6 7])),'all')
    error('A required original column contains NaN/Inf. Fix data provenance; never silently delete middle rows.');
end
fs=cfg.fs; f0=cfg.rollHz; N=size(A,1); t=(0:N-1)'/fs;
if N < 2000, error('At least 2 seconds are needed for the planned window audit.'); end
[~,stem]=fileparts(cfg.file);
out=fullfile(cfg.outputRoot,stem,'E04_REAL_REFERENCE');
if ~isfolder(out),mkdir(out);end
fprintf('E04: N=%d, assumed fs=%g Hz, reference fit f0=%g Hz\n',N,fs,f0);

%% Global descriptive statistics: evidence, not physical interpretation
names={'outer_frame_rate';'middle_frame_rate';'inner_frame_rate';...
       'gyro_a_raw';'column5_unverified';'accelerometer_col6';'gyro_s_raw'};
S=nan(7,10);
for c=1:7
    v=A(:,c);
    if any(~isfinite(v)), continue; end
    [amp,ph,r2]=localTone(t,v,f0);
    [amp6,~,r26]=localTone(t,v,2*f0);
    [fpeak,peakMag]=localPeak(v,fs,40);
    S(c,:)=[mean(v),std(v),median(v),min(v),max(v),amp,ph,r2,r26,fpeak];
end
T=array2table(S,'VariableNames',{'mean','std','median','minimum','maximum',...
                 'amplitude_f0','phase_f0_deg','fit_R2_f0','fit_R2_2f0','dominant_frequency_0_40Hz'});
T=[table((1:7)','VariableNames',{'column'}),table(string(names),'VariableNames',{'reported_role'}),T];
writetable(T,fullfile(out,'E04_column_statistics.csv'));

%% Independent roll-rate candidate SCREEN ONLY; not kinematic identification
% Since the expected roll frequency is ~3 Hz, equivalent speed is 1080 deg/s.
% A gimbal encoder rate and the actual body roll can differ by geometry.
F=table((1:3)',abs(S(1:3,1)),S(1:3,2),...
    abs(abs(S(1:3,1))-360*f0)/(360*f0),...
    S(1:3,2)./max(abs(S(1:3,1)),eps),...
    'VariableNames',{'column','absolute_mean_rate_reported_deg_s',...
    'std_rate_reported_deg_s','distance_fraction_to_3Hz_roll','relative_rate_variation'});
F.candidate_screen_only = F.distance_fraction_to_3Hz_roll < 0.15 & F.relative_rate_variation < 0.10;
writetable(F,fullfile(out,'E04_rollrate_candidates.csv'));

%% Windowed carrier phase relationship between independent col6 and a/s
windowN=round(2.0*fs); stepN=round(0.5*fs);
starts=1:stepN:(N-windowN+1); nw=numel(starts);
W=nan(nw,14);
for k=1:nw
    ii=starts(k):(starts(k)+windowN-1); tt=t(ii);
    [ampA,phaseA,r2A]=localTone(tt,A(ii,4),f0);
    [ampS,phaseS,r2S]=localTone(tt,A(ii,7),f0);
    [amp6,phase6,r26]=localTone(tt,A(ii,6),f0);
    [~,~,r2six2f0]=localTone(tt,A(ii,6),2*f0);
    da=localWrap(phase6-phaseA); ds=localWrap(phase6-phaseS);
    if r26<0.35, da=NaN; ds=NaN; end
    W(k,:)=[tt(1),tt(end),ampA,ampS,amp6,phaseA,phaseS,phase6,...
        localWrap(phaseS-phaseA),da,ds,r2A,r2S,r26];
end
P=array2table(W,'VariableNames',{'start_s','end_s','a_carrier_amp','s_carrier_amp',...
    'col6_carrier_amp','a_phase_deg','s_phase_deg','col6_phase_deg',...
    's_minus_a_phase_deg','col6_minus_a_deg_if_fit','col6_minus_s_deg_if_fit',...
    'a_fit_R2','s_fit_R2','col6_fit_R2'});
P.col6_fit_R2_at_2f0=nan(nw,1);
for k=1:nw
    ii=starts(k):(starts(k)+windowN-1);
    [~,~,P.col6_fit_R2_at_2f0(k)]=localTone(t(ii),A(ii,6),2*f0);
end
writetable(P,fullfile(out,'E04_phase_windows.csv'));

%% Diagnostic plot only; physical axes intentionally not re-labeled
fig=figure('Visible','off','Position',[100 100 1200 800]);
subplot(3,1,1);
plot(t,A(:,4),'Color',[.14 .45 .78]);hold on;
plot(t,A(:,7),'Color',[.90 .38 .18]);grid on;
legend('a: col4','s: col7','Location','best');xlabel('s');ylabel('original sensor units');
title('Unfiltered orthogonal gyro channels');
subplot(3,1,2);plot(t,A(:,6),'k');grid on;xlabel('s');ylabel('original col6 units');
title('Accelerometer col6: unverified sensor axis/units');
subplot(3,1,3);hold on;
for c=1:3,plot(t,A(:,c),'DisplayName',sprintf('gimbal col%d',c));end
xlabel('s');ylabel('reported deg/s');title('Turntable outer/middle/inner frame rates');
grid on;legend('Location','best');
saveas(fig,fullfile(out,'E04_original_columns.png'));close(fig);

%% Explicit limitations in report; do not automatically infer roll angle
fid=fopen(fullfile(out,'E04_interpretation_notes.txt'),'w');
fprintf(fid,'E04 real reference audit, no filtering / no parameter fitting.\n');
fprintf(fid,'Source: %s\nN=%d fs_assumed=%.6f Hz f0_assumed=%.6f Hz\n',cfg.file,N,fs,f0);
fprintf(fid,'Gimbal rate columns 1-3 may NOT directly equal body roll-rate.\n');
fprintf(fid,'Col6 is ONE measured acceleration channel; correlation at 3Hz is NOT proof of absolute roll phase.\n');
fprintf(fid,'No roll phase estimate, signed Omega truth or Euler-angle truth is constructed here.\n');
fprintf(fid,'Physical rotation-axis definitions, units, zero alignment and sign conventions are needed before reconstructing the reference attitude.\n');
fprintf(fid,'Observed 3 Hz fit R2: gyro a %.4f, gyro s %.4f, accel col6 %.4f; accel 6 Hz fit R2 %.4f.\n',...
    S(4,8),S(7,8),S(6,8),S(6,9));
fprintf(fid,'Gimbal rate-screen candidates (not proven axes): %s\n',mat2str(find(F.candidate_screen_only)'));
fprintf(fid,'Do NOT choose any method parameter using an already-viewed validation interval.\n');
fclose(fid);
fprintf('E04 complete: %s\n',out);

function [A,phaseDeg,R2] = localTone(t,x,f)
    t=double(t(:));x=double(x(:));
    C=[ones(numel(t),1),cos(2*pi*f*t),sin(2*pi*f*t)];
    b=C\x;
    A=hypot(b(2),b(3));
    phaseDeg=atan2d(-b(3),b(2));
    ss=sum((x-mean(x)).^2);
    if ss < 1e-20, R2=NaN; else,R2=1-sum((x-C*b).^2)/ss;end
end
function [fpeak,peak]=localPeak(x,fs,fmax)
    x=x(:)-mean(x); n=numel(x); X=abs(fft(x));
    inds=2:(min(floor(fmax*n/fs)+1,floor(n/2)+1));
    if isempty(inds),fpeak=NaN;peak=NaN;return;end
    [peak,k]=max(X(inds)); fpeak=(inds(k)-1)*fs/n;
end
function a=localWrap(a)
    a=mod(a+180,360)-180;
end
