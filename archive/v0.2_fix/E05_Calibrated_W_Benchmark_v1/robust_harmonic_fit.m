function fit = robust_harmonic_fit(angleDeg, y, cfg)
%ROBUST_HARMONIC_FIT 拟合 y=b+C*cos(phi)+S*sin(phi)。
% 截距置信区间基于逐滚转周期截距，避免把强自相关采样点当独立样本。

valid = isfinite(angleDeg) & isfinite(y);
angleDeg = angleDeg(valid);
y = y(valid);
theta = deg2rad(angleDeg);
X = [ones(size(theta)), cos(theta), sin(theta)];
base = robust_irls(X, y, cfg);

unwrapped = rad2deg(unwrap(theta));
progress = abs(unwrapped - unwrapped(1)) / 360;
cycleIndex = floor(progress);
cycleValues = unique(cycleIndex);
cycleBias = nan(numel(cycleValues), 1);
kept = false(numel(cycleValues), 1);
for k = 1:numel(cycleValues)
    mask = cycleIndex == cycleValues(k);
    phaseSpan = max(unwrapped(mask)) - min(unwrapped(mask));
    if sum(mask) >= 12 && phaseSpan >= 250
        oneCycle = robust_irls(X(mask, :), y(mask), cfg);
        cycleBias(k) = oneCycle.beta(1);
        kept(k) = true;
    end
end
cycleBias = cycleBias(kept);
if numel(cycleBias) >= 4
    [biasLocation, ciLow, ciHigh, biasSigma] = robust_location_ci(cycleBias, cfg);
else
    biasLocation = base.beta(1);
    standardError = sqrt(max(base.covariance(1, 1), 0));
    ciLow = biasLocation - cfg.confidenceZ * standardError;
    ciHigh = biasLocation + cfg.confidenceZ * standardError;
    biasSigma = base.robustScale;
end

fit.bias = biasLocation;
fit.biasCiLow = ciLow;
fit.biasCiHigh = ciHigh;
fit.biasSigma = biasSigma;
fit.cosCoefficient = base.beta(2);
fit.sinCoefficient = base.beta(3);
fit.complexCoefficient = complex(base.beta(2), base.beta(3));
fit.amplitude = abs(fit.complexCoefficient);
fit.phaseDeg = atan2d(fit.sinCoefficient, fit.cosCoefficient);
fit.r2 = base.r2;
fit.residualScale = base.robustScale;
fit.outlierFraction = base.outlierFraction;
fit.cycleCount = numel(cycleBias);
fit.prediction = X * [fit.bias; fit.cosCoefficient; fit.sinCoefficient];
fit.residual = y - fit.prediction;
end

function [location, ciLow, ciHigh, spread] = robust_location_ci(values, cfg)
values = values(isfinite(values));
location = median(values);
for iteration = 1:cfg.huberMaxIterations
    spread = 1.4826 * median(abs(values - location));
    if spread <= eps
        spread = std(values);
    end
    if spread <= eps
        break;
    end
    cutoff = cfg.huberTune * spread;
    weights = min(1, cutoff ./ max(abs(values - location), eps));
    updated = sum(weights .* values) / sum(weights);
    if abs(updated - location) <= cfg.huberConvergenceTolerance * (1 + abs(location))
        location = updated;
        break;
    end
    location = updated;
end
if ~exist('spread', 'var') || ~isfinite(spread)
    spread = 0;
end
standardError = std(values) / sqrt(max(numel(values), 1));
ciLow = location - cfg.confidenceZ * standardError;
ciHigh = location + cfg.confidenceZ * standardError;
end

