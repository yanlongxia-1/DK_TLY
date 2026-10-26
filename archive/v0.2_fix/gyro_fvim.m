function [bestK, bestAlpha, history, details] = gyro_fvim(xTrain,cfg,strategy)
% Four-leader FVIM-style optimizer with documented reproducible choices.
% Paper leaves initialization/threshold/fitness details incomplete.
% This is NOT claimed to reproduce author source code exactly.
rng(cfg.seed,'twister');
popN=max(6,cfg.fvim.population); nIter=cfg.fvim.iterations;
low=[cfg.searchK(1),log10(cfg.searchAlpha(1))];
high=[cfg.searchK(2),log10(cfg.searchAlpha(2))];
pop=low+rand(popN,2).*(high-low);
% Deterministic coverage of extremes and center before random exploration.
pop(1,:)=low; pop(2,:)=high;
pop(3,:)=(low+high)/2; pop(4,:)=[low(1),high(2)];
scores=inf(popN,1);
history=zeros(nIter+1,4); calls=0;
for i=1:popN
    scores(i)=objective(pop(i,:));
end
[score,ix]=min(scores);
history(1,:)=[0,round(pop(ix,1)),10^pop(ix,2),score];
stall=0; prevBest=min(scores);
for it=1:nIter
    [~,idx]=sort(scores);
    leaders=pop(idx(1:4),:);
    amp=2*(1-it/(nIter+1));
    for i=1:popN
        moves=zeros(4,2);
        for j=1:4
            A=amp*(2*rand(1,2)-1); C=2*rand(1,2);
            dist=abs(C.*leaders(j,:)-pop(i,:));
            moves(j,:)=leaders(j,:)-A.*dist;
        end
        trial=mean(moves,1);
        if i>ceil(popN/2) && rand<0.35
            trial=trial+0.12*randn(1,2).*(high-low);
        end
        if stall>=3 && i>floor(popN/2)
            % Restart part of the population if global-best fitness stalls.
            trial=low+rand(1,2).*(high-low);
        end
        trial=max(low,min(high,trial));
        newScore=objective(trial);
        if newScore<scores(i)
            pop(i,:)=trial; scores(i)=newScore;
        end
    end
    [score,ix]=min(scores);
    if score<prevBest-1e-10
        stall=0;prevBest=score;
    else
        stall=stall+1;
    end
    history(it+1,:)=[it,round(pop(ix,1)),10^pop(ix,2),score];
end
[~,ix]=min(scores);
bestK=round(pop(ix,1));bestAlpha=10^pop(ix,2);
details.calls=calls;
details.objective='Min envelope entropy (operational surrogate); physical mode adds amplitude/phase/band preservation penalties';
details.bestFitness=min(scores);
details.boundaryK=(bestK<=cfg.searchK(1) || bestK>=cfg.searchK(2));
details.bestK=bestK;details.bestAlpha=bestAlpha;
if ~isfinite(details.bestFitness) || details.bestFitness>=1e6
    error('gyro_fvim:noValidCandidate','No converged MVMD candidate found on TRAIN data; do not evaluate validation.');
end

    function score=objective(v)
        calls=calls+1;
        K=max(cfg.searchK(1),min(cfg.searchK(2),round(v(1))));
        alpha=10^v(2);
        try
            [modes,centers,mvinfo]=gyro_mvmd(xTrain,cfg.fs,K,alpha, ...
                cfg.mvmd.searchIterations,cfg.mvmd.tolerance,cfg.mvmd.tau);
            if ~mvinfo.converged
                score=1e7+1000*mvinfo.delta; return;
            end
            h=zeros(K,2);
            for kk=1:K
                for cc=1:2
                    env=abs(gyro_analytic(modes(:,cc,kk)));
                    p=env/(sum(env)+eps);
                    h(kk,cc)=-sum(p.*log(p+eps))/log(numel(p));
                end
            end
            % Documented proxy for paper's incompletely specified min envelope entropy.
            score=mean(min(h,[],1));
            % Avoid degenerately empty modes that spuriously yield low envelope entropy.
            modeEnergy=squeeze(sum(sum(modes.^2,1),2));
            energyShare=modeEnergy/(sum(modeEnergy)+eps);
            if any(energyShare<1e-5),score=score+1;end
            if strcmpi(strategy,'physics')
                [r,~]=gyro_reconstruct(modes,centers,xTrain,cfg,'physics');
                u=gyro_sinefit(xTrain,cfg.fs,cfg.rollHz);
                z=gyro_sinefit(r,cfg.fs,cfg.rollHz);
                aerr=max(abs(z.amplitude./(u.amplitude+eps)-1));
                ur=angle(exp(1i*(u.phase(2)-u.phase(1))));
                zr=angle(exp(1i*(z.phase(2)-z.phase(1))));
                pd=abs(rad2deg(angle(exp(1i*(zr-ur)))));
                er=band_power_loss(xTrain,r,cfg.fs,cfg.rollHz,cfg.physics.bandHz);
                w=cfg.physics.weights;
                score=score+w(1)*max(0,(aerr-cfg.physics.ampTol)/max(cfg.physics.ampTol,eps))^2 ...
                    +w(2)*max(0,(pd-cfg.physics.phaseTolDeg)/max(cfg.physics.phaseTolDeg,eps))^2 ...
                    +w(3)*max(0,er/0.10)^2;
            end
            if ~isfinite(score),score=1e9;end
        catch
            score=1e9;
        end
    end
end

function loss=band_power_loss(x,z,fs,f0,halfWidth)
N=size(x,1);
f=(0:N-1)'*fs/N;
mask=(f>=max(0,f0-halfWidth)&f<=f0+halfWidth) | ...
    (f>=fs-(f0+halfWidth)&f<=fs-max(0,f0-halfWidth));
% Energy loss estimate in the signal vicinity; raw is a noisy proxy.
A=fft(x-mean(x,1),[],1); B=fft(z-mean(z,1),[],1);
ea=sum(abs(A(mask,:)).^2,'all');eb=sum(abs(B(mask,:)).^2,'all');
loss=max(0,(ea-eb)/(ea+eps));
end
