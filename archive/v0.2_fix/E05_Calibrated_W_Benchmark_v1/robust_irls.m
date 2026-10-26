function fit = robust_irls(X, y, cfg)
%ROBUST_IRLS 无工具箱依赖的 Huber 迭代加权最小二乘。

valid = all(isfinite(X), 2) & isfinite(y);
X = X(valid, :);
y = y(valid);
if size(X, 1) <= size(X, 2)
    error('robust_irls:InsufficientData', '稳健回归的有效样本不足。');
end

beta = X \ y;
weights = ones(size(y));
for iteration = 1:cfg.huberMaxIterations
    residual = y - X * beta;
    center = median(residual);
    scale = 1.4826 * median(abs(residual - center));
    if ~isfinite(scale) || scale <= eps(max(abs(y)) + 1)
        scale = sqrt(mean(residual .^ 2) + eps);
    end
    cutoff = cfg.huberTune * scale;
    weights = min(1, cutoff ./ max(abs(residual - center), eps));
    sqrtWeights = sqrt(weights);
    betaNew = (X .* sqrtWeights) \ (y .* sqrtWeights);
    if norm(betaNew - beta) <= cfg.huberConvergenceTolerance * (1 + norm(beta))
        beta = betaNew;
        break;
    end
    beta = betaNew;
end

residual = y - X * beta;
robustScale = 1.4826 * median(abs(residual - median(residual)));
if robustScale <= eps
    robustScale = sqrt(mean(residual .^ 2) + eps);
end
weightedSse = sum(weights .* residual .^ 2);
degreesFreedom = max(1, size(X, 1) - size(X, 2));
covariance = (weightedSse / degreesFreedom) * pinv(X' * (X .* weights));
sst = sum((y - mean(y)) .^ 2);

fit.beta = beta;
fit.residual = residual;
fit.weights = weights;
fit.robustScale = robustScale;
fit.covariance = covariance;
fit.r2 = 1 - sum(residual .^ 2) / max(sst, eps);
fit.outlierFraction = mean(abs(residual - median(residual)) > ...
    cfg.outlierMadThreshold * robustScale);
fit.iterations = iteration;
end

