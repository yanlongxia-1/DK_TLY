function [yhat, details] = gyro_reconstruct(modes, centerHz, x, cfg, strategy)
% Classification uses BOTH channels jointly, and keeps mode pairing.
% IMPORTANT: 'mixed' requires BOTH sufficient cosine similarity AND not-too-high entropy.
% Entropy thresholds are preregistered engineering settings, NOT exact paper source code.
% cs/SampEn thresholds are explicit operational assumptions; not author code.
[N,C,K]=size(modes);
if C~=2, error('This lab expects a/s dual-channel data'); end
classes=zeros(K,1); % 0=discard; 1=mixed-wavelet; 2=keep
cs=zeros(K,1); sp=zeros(K,1); carrier=zeros(K,1);
bandShares=nan(K,1); modeEnergy=zeros(K,1);
t=(0:N-1)'/cfg.fs;
phasor=exp(-1i*2*pi*cfg.rollHz*t);
carrierSignal=max(abs(sum((x-mean(x,1)).*phasor,1)),eps);
for k=1:K
    m=modes(:,:,k);
    modeEnergy(k)=sum(m(:).^2);
    dot12=sum(m.*(x-mean(x,1)),1);
    den=sqrt(sum(m.^2,1).*sum((x-mean(x,1)).^2,1))+eps;
    cs(k)=mean(abs(dot12./den));
    sj=zeros(1,2);
    for c=1:2
        sj(c)=gyro_sampen(m(:,c),cfg.class.maxSampPoints,cfg.class.sampTolerance);
    end
    sp(k)=mean(sj);
    carrier(k)=mean(abs(sum(m.*phasor,1))./carrierSignal);
    if cs(k)>=cfg.class.csKeep && sp(k)<=cfg.class.sampKeep
        classes(k)=2;
    elseif cs(k)>=cfg.class.csMixed && sp(k)<=cfg.class.sampMixed
        classes(k)=1;
    end
end
beforePhysics=classes;
if strcmpi(strategy,'physics')
    % Paired band-energy protection: protect plausible carrier/sideband energy
    % only when BOTH channels carry energy in the predeclared physical band.
    bandShares=zeros(K,1);
    f=(0:N-1)'*cfg.fs/N;
    mask=(f>=max(0,cfg.rollHz-cfg.physics.bandHz) & f<=cfg.rollHz+cfg.physics.bandHz) | ...
         (f>=cfg.fs-(cfg.rollHz+cfg.physics.bandHz) & f<=cfg.fs-max(0,cfg.rollHz-cfg.physics.bandHz));
    base=abs(fft(x-mean(x,1),[],1)).^2;
    basePow=sum(base(mask,:),1)+eps;
    for k=1:K
        mm=abs(fft(modes(:,:,k),[],1)).^2;
        bandShares(k)=min(sum(mm(mask,:),1)./basePow);
    end
    for k=1:K
        if abs(centerHz(k)-cfg.rollHz)<=cfg.physics.protectHz && ...
                carrier(k)>=cfg.physics.minCarrierShare
            classes(k)=2; % guard both channels, not independently
        elseif centerHz(k)<=cfg.rollHz+cfg.physics.bandHz && ...
                bandShares(k)>=cfg.physics.minBandShare && classes(k)<2
            classes(k)=2; % possible motion-modulation sideband; verified jointly
        end
    end
elseif ~strcmpi(strategy,'paper')
    error('Unknown strategy: %s',strategy);
end
if ~any(classes)
    [~,k]=max(cs+0.5*carrier);
    classes(k)=2;
end
yhat=zeros(N,C);
for k=1:K
    if classes(k)==2
        yhat=yhat+modes(:,:,k);
    elseif classes(k)==1
        for c=1:C
            yhat(:,c)=yhat(:,c)+gyro_haar(modes(:,c,k),cfg.wavelet.levels,...
                cfg.wavelet.thresholdScale);
        end
    end
end
% DC is preserved, not discarded. Manual calibration is done before this call.
yhat=yhat+mean(x,1);
details=table((1:K)',centerHz(:),cs,sp,carrier, ...
    modeEnergy/(sum(modeEnergy)+eps),beforePhysics,bandShares, ...
    double(classes>beforePhysics),classes, ...
    'VariableNames',{'mode','center_Hz','cosine_abs','sampEn','carrier_share', ...
    'energy_share','classification_before_physics','physics_band_share', ...
    'physical_protection_applied','class_0drop_1mix_2keep'});
end
