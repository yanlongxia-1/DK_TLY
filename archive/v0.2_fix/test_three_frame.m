clear; clc; close all;

% 修改为你的真实文件路径
data = 'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验\20200824-1-2\W3-20.txt';

fs = 1000;
N = size(data, 1);
t = (0:N-1)' / fs;

% 假设第1、3列为角位置（单位：deg）
outerRaw = data(:,1);
innerRaw = data(:,3);

% 消除0/360度及-360/0度回绕
outerAngle = rad2deg(unwrap(deg2rad(outerRaw)));
innerAngle = rad2deg(unwrap(deg2rad(innerRaw)));

% 计算回绕次数
outerWraps = sum(abs(diff(outerRaw)) > 180);
innerWraps = sum(abs(diff(innerRaw)) > 180);

% 线性拟合得到平均角速度
p1 = polyfit(t, outerAngle, 1);
p3 = polyfit(t, innerAngle, 1);

% 拟合残差与拟合优度
fitOuter = polyval(p1,t);
fitInner = polyval(p3,t);

R2_outer = 1 - sum((outerAngle-fitOuter).^2) / ...
                sum((outerAngle-mean(outerAngle)).^2);

R2_inner = 1 - sum((innerAngle-fitInner).^2) / ...
                sum((innerAngle-mean(innerAngle)).^2);

fprintf('\n===== W3-20 角位置验证 =====\n');
fprintf('外框回绕次数: %d\n', outerWraps);
fprintf('内框回绕次数: %d\n', innerWraps);
fprintf('外框拟合角速度: %.6f deg/s\n', p1(1));
fprintf('内框拟合角速度: %.6f deg/s\n', p3(1));
fprintf('内框滚转频率绝对值: %.6f Hz\n', abs(p3(1))/360);
fprintf('外框线性拟合R2: %.8f\n', R2_outer);
fprintf('内框线性拟合R2: %.8f\n', R2_inner);

figure;
subplot(2,1,1);
plot(t,outerAngle,'b',t,fitOuter,'r--');
title('外框角度展开与线性拟合');
xlabel('Time(s)'); ylabel('Angle(deg)');
grid on;

subplot(2,1,2);
plot(t,innerAngle,'b',t,fitInner,'r--');
title('内框角度展开与线性拟合');
xlabel('Time(s)'); ylabel('Angle(deg)');
grid on;