function E07D_selftest(cfg)
%E07D_SELFTEST 无随机数据与完整FVIM寻优的快速确定性检验。
cs=[.4;.321768;.321331;.01];
sp=[.8;1.793180;1.815064;2.3];
[k,m]=E07D_soft_weights(cs,sp,cfg);
assert(all(isfinite(k+m)) && all(k>=0 & k<=1) && all(m>=0 & m<=1));
assert(k(1)>.5,'强有效模态应有较大保留权重');
assert(abs(m(2)-m(3))<.15,'阈值附近应连续变化而非硬切换');
assert(m(4)<.01,'高熵且低CS噪声模态不应保留');
t=(0:999)'/1000;
x=[cos(2*pi*3*t),sin(2*pi*3*t)];
assert(isequal(size(x),[1000 2]));
fprintf('E07-D SELFTEST PASS: 连续权重及基本数组尺寸检查通过。\n');
end
