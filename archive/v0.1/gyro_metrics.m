function T = gyro_metrics(raw, filt, fs, f0, name)
% Metrics relative to MEASUREMENTS, never mislabeled ground-truth errors.
a=gyro_sinefit(raw,fs,f0);
b=gyro_sinefit(filt,fs,f0);
phaseShift=rad2deg(angle(exp(1i*(b.phase-a.phase))));
pa=rad2deg(angle(exp(1i*(a.phase(2)-a.phase(1)))));
pb=rad2deg(angle(exp(1i*(b.phase(2)-b.phase(1)))));
relativeChange=rad2deg(angle(exp(1i*deg2rad(pb-pa))));
N=size(raw,1);
chan={'a_col4';'s_col7'};
method=repmat({name},2,1);
RMSE=zeros(2,1); NRMSE=zeros(2,1); R=zeros(2,1); snrGain=zeros(2,1);
meanOffset=zeros(2,1);
roughnessKeep=zeros(2,1);residualPowerRatio=zeros(2,1);noncarrierDrop=zeros(2,1);
for c=1:2
    RMSE(c)=sqrt(mean((filt(:,c)-raw(:,c)).^2));
    NRMSE(c)=100*RMSE(c)/(std(raw(:,c))+eps);
    corrmat=corrcoef(raw(:,c),filt(:,c));
    if numel(corrmat)==4, R(c)=corrmat(1,2); else, R(c)=NaN; end
    snrGain(c)=b.snrProxy(c)-a.snrProxy(c);
    meanOffset(c)=mean(filt(:,c)-raw(:,c));
    roughnessKeep(c)=100*sum(diff(filt(:,c),2).^2)/(sum(diff(raw(:,c),2).^2)+eps);
    residualPowerRatio(c)=10*log10((sum((filt(:,c)-mean(filt(:,c))).^2)+eps)/...
        (sum((raw(:,c)-filt(:,c)).^2)+eps));
    noncarrierDrop(c)=10*log10((sum(a.residual(:,c).^2)+eps)/(sum(b.residual(:,c).^2)+eps));
end
T=table(method,chan,RMSE,NRMSE,R,100*b.amplitude(:)./(a.amplitude(:)+eps), ...
    phaseShift(:),phaseShift(:)/(360*f0)*1000, a.snrProxy(:),b.snrProxy(:),snrGain,...
    meanOffset,roughnessKeep,residualPowerRatio,noncarrierDrop,repmat(pa,2,1),repmat(pb,2,1),repmat(relativeChange,2,1), ...
    repmat(N,2,1), ...
    'VariableNames',{'method','channel','RMSE_vs_measured','NRMSE_vs_std_pct', ...
    'correlation','carrier_amp_keep_pct','phase_shift_deg','equiv_delay_ms', ...
    'raw_snr_proxy_dB','filter_snr_proxy_dB','snr_proxy_gain_dB', ...
    'mean_offset','second_diff_roughness_keep_pct','output_to_residual_power_dB',...
    'noncarrier_residual_drop_dB','raw_relative_phase_deg','filt_relative_phase_deg', ...
    'relative_phase_change_deg','nSamples'});
end
