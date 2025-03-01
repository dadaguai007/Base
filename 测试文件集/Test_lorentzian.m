clc;clear;close;
% 定义洛伦兹函数
lorentzian = @(x, a, b, c, d) a + b ./ ((x - c).^2 + d^2);

% 准备数据
xData = linspace(-10, 10, 100);  % 生成X数据
a = 10; b = 10; c = 0; d = 1;  % 洛伦兹函数参数
yData = lorentzian(xData, a, b, c, d);  % 生成Y数据

% 绘制洛伦兹曲线
figure;
plot(xData, yData, 'LineWidth', 2);
title('洛伦兹曲线');
xlabel('X轴');
ylabel('Y轴');
grid on;



