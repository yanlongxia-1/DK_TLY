function [metrics, seriesOut] = E09A0_calculate(encoderAngleAllDeg, rates, cfg)
%E09A0_CALCULATE 使用编码器角度真值的外框积分角度对照。
% 重要：不是机体pitch/yaw角度，任何方法都不能用评价区间的编码器角度修正零偏/滤波参数。
assert(isvector(encoderAngleAllDeg) && isnumeric(encoderAngleAllDeg));
assert(istable(rates) && ismember('time_s',rates.Properties.VariableNames));
assert(ismember('reference_outer_deg_s',rates.Properties.VariableNames), ...
    '找不到E08-K参考转速列，必须读取正确的W类序列CSV');
t = double(rates.time_s(:));
assert(numel(t)>50 && all(isfinite(t)) && all(diff(t)>0),'时间戳无效');
assert(all(abs(diff(t)-1/cfg.fs)<cfg.sampleTimeToleranceSec),...
    'CSV时间戳不是连续的1/fs采样：请检查丢样和时间基');
idxFloat = t*cfg.fs+1;
ix = round(idxFloat);
assert(all(abs(idxFloat-ix)/cfg.fs < cfg.sampleTimeToleranceSec),...
    'CSV时间戳与原始TXT样本索引不一致');
assert(all(ix>=1) && all(ix<=numel(encoderAngleAllDeg)),...
    'CSV引用的样本超出TXT文件范围');
rawEncoder = double(encoderAngleAllDeg(:));
assert(all(isfinite(rawEncoder)),'转台角度存在NaN/Inf');
encoderUnwrapped = rad2deg(unwrap(deg2rad(rawEncoder)));
ref = encoderUnwrapped(ix);
% 同一评价起点；只使用第1个编码器样本确定初始角度，不使用后续真值拟合误差。
vfit = double(rates.reference_outer_deg_s(:));
assert(all(isfinite(vfit)),'参考角速度列存在NaN/Inf');
refFit = ref(1)+cumtrapz(t,vfit);
seriesOut=table(t,ix,ref,refFit,'VariableNames',...
    {'time_s','originalSampleIndex','encoder_outer_angle_deg','fitted_outer_angle_deg'});
metrics=table();
reserved={'time_s','reference_outer_deg_s'};
vars=rates.Properties.VariableNames;
for i=1:numel(vars)
    name=vars{i};
    if any(strcmp(name,reserved)),continue;end
    v = double(rates.(name));
    if numel(v)~=numel(t) || any(~isfinite(v)),continue;end
    estimate = ref(1) + cumtrapz(t, v(:));
    errorEncoder = estimate-ref;
    errorFitted = estimate-refFit;
    endpoint = errorEncoder(end);
    Trow=table(string(name),numel(t),t(end)-t(1), ...
        sqrt(mean(errorEncoder.^2)),mean(abs(errorEncoder)),max(abs(errorEncoder)), ...
        endpoint, sqrt(mean(errorFitted.^2)),errorFitted(end),mean(v-vfit), ...
        'VariableNames', {'method','samples','duration_s', ...
        'RMSE_encoder_deg','MAE_encoder_deg','MaxAbsError_encoder_deg', ...
        'EndError_encoder_deg','RMSE_fitted_deg','EndError_fitted_deg','RateBias_vs_fitted_deg_s'});
    if isempty(metrics),metrics=Trow;else,metrics=[metrics;Trow];end %#ok<AGROW>
    seriesOut.(matlab.lang.makeValidName(['angle_' name]))=estimate;
    seriesOut.(matlab.lang.makeValidName(['error_' name]))=errorEncoder;
end
assert(~isempty(metrics),'没有任何有效的角速度方法列');
end
