%% CAL01: fast, train-only harmonic effective-offset/gain calibration candidates
% Input TXT col1/2/3 = OUTER/MIDDLE/INNER angles (deg).
% col4=a; col7=s. Output parameters have the ORIGINAL data units.
% Estimates are NOT certified zero-rate electronic bias or true sensitivity
% until calibration geometry and stationary (zero-input) checks are established.
clear; clc; close all;
cfg = CAL01_config();
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
header = {'case_name','channel','N','train_n','check_n','f_roll_set_Hz', ...
 'outer_nominal_deg_s','outer_measured_deg_s','outer_angle_linear_R2', ...
 'inner_measured_Hz_signed','middle_abs_median_deg','bias_effective_train', ...
 'carrier_amplitude_train','K_eff_using_outer_rate', ...
 'K_nominal_using_filename','carrier_phase_train_deg', ...
 'train_harmonic_R2','train_residual_RMS', ...
 'bias_check_only','carrier_amplitude_check_only', ...
 'K_eff_check_only','relative_bias_drift','relative_gain_drift_pct', ...
 'phase_s_minus_a_train_deg','geometry_verified'};
rows = cell(0,numel(header));
notes = {};
for ic = 1:numel(cfg.cases)
 cs = cfg.cases(ic);
 if ~isfile(cs.file)
   warning('Missing %s: %s',cs.name,cs.file);
   notes{end+1,1}=sprintf('MISSING: %s (%s)',cs.name,cs.file); %#ok<SAGROW>
   continue
 end
 M=readmatrix(cs.file,'FileType','text');
 if size(M,2)<7,error('%s needs >=7 numeric columns',cs.name);end
 D=double(M(:,[1,2,3,4,7]));
 if any(~isfinite(D(:))),error('%s contains NaN/Inf in columns 1/2/3/4/7',cs.name);end
 N=size(M,1); t=(0:N-1)'/cfg.fs;
 train=segment(N,cfg.trainFraction); check=segment(N,cfg.checkFraction);
 if isempty(train)||isempty(check)||~isempty(intersect(train,check)),error('Invalid train/check split');end
 outer=rad2deg(unwrap(deg2rad(double(M(:,1)))));
 inner=rad2deg(unwrap(deg2rad(double(M(:,3)))));
 middle=double(M(:,2));
 pOuter=polyfit(t(train),outer(train),1);
 outerRate=abs(pOuter(1));
 outerR2=oneMinusSSE(outer(train),polyval(pOuter,t(train)));
 pInner=polyfit(t(train),inner(train),1);
 innerFreq=pInner(1)/360;
 middleMed=median(abs(middle(train)));
 if outerRate<1e-7,error('%s: outer rate zero; cannot estimate gain',cs.name);end
 y=double(M(:,[4,7]));
 fitTr=cell(1,2);fitChk=cell(1,2);
 for k=1:2
   fitTr{k}=harmonicFit(t(train),y(train,k),cs.fRollHz);
   fitChk{k}=harmonicFit(t(check),y(check,k),cs.fRollHz);
 end
 phaseDiff=wrap180(fitTr{2}.phaseDeg-fitTr{1}.phaseDeg);
 for k=1:2
  F=fitTr{k};V=fitChk{k};
  Ktrain=F.amplitude/outerRate;
  Kcheck=V.amplitude/outerRate; % diagnostic ONLY; locked outer slope from TRAIN
  biasDrift=V.bias-F.bias;
  gainDriftPct=100*(Kcheck/Ktrain-1);
  chanNames={'a','s'};
  rows(end+1,:) = {cs.name, chanNames{k}, N, numel(train), numel(check), cs.fRollHz, ...
    cs.nominalOuterDegS,outerRate,outerR2,innerFreq,middleMed,F.bias, ...
    F.amplitude,Ktrain,F.amplitude/cs.nominalOuterDegS,F.phaseDeg, ...
    F.R2,F.rms,V.bias,V.amplitude,Kcheck,biasDrift,gainDriftPct, ...
    phaseDiff, double(cfg.geometryVerified)}; %#ok<SAGROW>
 end
 if ~cfg.geometryVerified
   notes{end+1,1}=sprintf('%s: K_eff is amplitude/outer-angle-slope, NOT VERIFIED physical scale K.',cs.name); %#ok<SAGROW>
 end
 if abs(outerRate-cs.nominalOuterDegS)>max(0.7,0.06*cs.nominalOuterDegS)
   notes{end+1,1}=sprintf('%s: measured outer rate %.3f differs from nominal %.3f deg/s',cs.name,outerRate,cs.nominalOuterDegS); %#ok<SAGROW>
 end
 if outerR2<0.995
   notes{end+1,1}=sprintf('%s: outer angle slope R2=%.4f low; inspect angular reference',cs.name,outerR2); %#ok<SAGROW>
 end
 if abs(abs(innerFreq)-cs.fRollHz)>0.05
   notes{end+1,1}=sprintf('%s: inner frame frequency %.4fHz differs from set %.3fHz',cs.name,innerFreq,cs.fRollHz); %#ok<SAGROW>
 end
 if abs(abs(phaseDiff)-90)>5
   notes{end+1,1}=sprintf('%s: channel relative phase %.2fdeg not close to quadrature',cs.name,phaseDiff); %#ok<SAGROW>
 end
 if middleMed>2
   notes{end+1,1}=sprintf('%s: middle angle median absolute %.2fdeg; do not assume Omega_transverse=outer rate',cs.name,middleMed); %#ok<SAGROW>
 end
 fig=figure('Visible','off','Color','w');
 for k=1:2
   subplot(2,1,k);
   plot(t(train),y(train,k),'Color',[0.7,0.7,0.7]);hold on;
   plot(t(train),predictHarmonic(t(train),fitTr{k},cs.fRollHz),'b','LineWidth',1.3);
   ylabel(sprintf('Channel %s',char('a'+18*(k-1))));grid on;
   if k==1,title(sprintf('%s: TRAIN ONLY sinusoid + effective offset',cs.name),'Interpreter','none');end
   if k==2,xlabel('Time (s)');end
 end
 saveas(fig,fullfile(cfg.outputDir,[cs.name '_CAL01_fit.png']));close(fig);
 fprintf('\n%s: measured outer=%.5f deg/s, inner=%+.5fHz, s-a phase=%+.3f deg\n', ...
   cs.name,outerRate,innerFreq,phaseDiff);
 fprintf('  a: b_eff=%+.8f, A=%.8f, K_eff=%.9f\n',fitTr{1}.bias,fitTr{1}.amplitude,fitTr{1}.amplitude/outerRate);
 fprintf('  s: b_eff=%+.8f, A=%.8f, K_eff=%.9f\n',fitTr{2}.bias,fitTr{2}.amplitude,fitTr{2}.amplitude/outerRate);
end
if isempty(rows),error('No real TXT input loaded. Edit CAL01_config.m first.');end
T=cell2table(rows,'VariableNames',header);
for j=3:numel(header)
  col=T.(header{j});
  if all(cellfun(@(v) isnumeric(v)&&isscalar(v),col))
    T.(header{j})=cell2mat(col);
  end
end
writetable(T,fullfile(cfg.outputDir,'CAL01_parameters.csv'));
% Compare bias/gain consistency. These numbers are diagnostic, not proof
% that intercept is a static electrical bias across different conditions.
if numel(unique(T.case_name))>1
 channelNames={'a','s'};
 for k=1:2
  sub=T(strcmp(T.channel,channelNames{k}),:);
  if height(sub)>1
    fprintf('Across-file %s: effective bias range %.6f; K_eff range %.8f\n', ...
       sub.channel{1},max(sub.bias_effective_train)-min(sub.bias_effective_train), ...
       max(sub.K_eff_using_outer_rate)-min(sub.K_eff_using_outer_rate));
  end
 end
end
fid=fopen(fullfile(cfg.outputDir,'CAL01_E05_candidate_config.txt'),'w');
if fid>=0
 fprintf(fid,'CAL01 CANDIDATES ONLY. They are not independent laboratory calibrations.\n');
 fprintf(fid,'Do NOT use these to claim absolute attitude accuracy without confirmed gimbal geometry.\n\n');
 for ic=1:numel(cfg.cases)
  nm=cfg.cases(ic).name;
  a=T(strcmp(T.case_name,nm)&strcmp(T.channel,'a'),:);
  s=T(strcmp(T.case_name,nm)&strcmp(T.channel,'s'),:);
  if height(a)==1 && height(s)==1
    fprintf(fid,'%% %s - effective offset and apparent gain from TRAIN only\n',nm);
    fprintf(fid,'bias = [%.11g, %.11g];\n',a.bias_effective_train,s.bias_effective_train);
    fprintf(fid,'scale_candidate = [%.11g, %.11g]; %% NOT verified physical K\n\n',a.K_eff_using_outer_rate,s.K_eff_using_outer_rate);
  end
 end
 fclose(fid);
end
fid=fopen(fullfile(cfg.outputDir,'CAL01_NOTES.txt'),'w');
if fid>=0
 fprintf(fid,'TRAIN-only estimates. CHECK is diagnostic only; never refit on CHECK.\n');
 fprintf(fid,'Signals: col4=a and col7=s; gimbal col1-3 are ANGLES deg.\n');
 fprintf(fid,'Fit: x(t) = b_effective + P*sin(2*pi*f*t) + Q*cos(2*pi*f*t).\n');
 fprintf(fid,'K_eff = fitted carrier amplitude / abs(train outer angle slope).\n');
 fprintf(fid,'K_eff equals physical scale K only if transverse rate is established to equal outer gimbal rate.\n');
 fprintf(fid,'b_effective is a harmonic intercept, NOT verified static electronic bias.\n');
 fprintf(fid,'To certify b: acquire zero-transverse-input reference with matching roll/electronics conditions.\n');
 fprintf(fid,'For formal K: confirm gimbal-axis/mounting geometry and test >=2 known nonzero rates.\n');
 fprintf(fid,'No FVIM, MVMD or large parameter grid is run in CAL01.\n');
 for i=1:numel(notes),fprintf(fid,'WARNING: %s\n',notes{i});end
 fclose(fid);
end
fprintf('\nFinished. Output: %s\n',cfg.outputDir);

function ids=segment(N,f)
 ids=(max(1,round(f(1)*N)):min(N,round(f(2)*N)))';
end
function F=harmonicFit(t,y,f)
 t=t(:);y=y(:);
 X=[ones(numel(t),1),sin(2*pi*f*t),cos(2*pi*f*t)];
 beta=X\y;
 for iter=1:6
   r=y-X*beta; sigma=1.4826*median(abs(r-median(r)))+1e-12;
   w=min(1,(1.345*sigma)./max(abs(r),1e-12));
   W=sqrt(w); beta=(X.*W)\(y.*W);
 end
 residual=y-X*beta;
 F.bias=beta(1);F.amplitude=hypot(beta(2),beta(3));
 F.phaseDeg=atan2d(beta(3),beta(2));
 F.rms=sqrt(mean(residual.^2));
 F.R2=oneMinusSSE(y,X*beta);
 F.sinCoef=beta(2);F.cosCoef=beta(3);
end
function y=predictHarmonic(t,F,f)
 y=F.bias+F.sinCoef*sin(2*pi*f*t)+F.cosCoef*cos(2*pi*f*t);
end
function r2=oneMinusSSE(y,fit)
 denom=sum((y-mean(y)).^2);
 if denom<1e-20,r2=NaN;else,r2=1-sum((y-fit).^2)/denom;end
end
function v=wrap180(v)
 v=mod(v+180,360)-180;
end
