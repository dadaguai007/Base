% Test Hpuls_shape
clc;clear;
% 设置采样率
fs = 1000; % 采样率，单位为Hz

% 创建时间向量，从0到脉冲持续时间
t = 0:1/fs:3/fs;

% 创建一个矩形脉冲形状，持续时间为4个样本
h_rect = rectpuls(t,3/fs);

% 创建一个Raised Cosine脉冲形状，持续时间为4个样本，滚降因子为0.5
h_rc = rcosdesign(0.5, 4/fs, fs);

% 创建频率向量，从0到采样率的一半
f = 0:fs/2;

% 计算矩形脉冲的频率响应
H_rect = Hpulsshape(f, h_rect, fs);

% 计算Raised Cosine脉冲的频率响应
H_rc = Hpulsshape(f, h_rc, fs);

% 绘制频率响应的幅度
figure;
subplot(2,1,1);
plot(f, abs(H_rect));
title('矩形脉冲的频率响应幅度');
xlabel('频率 (Hz)');
ylabel('幅度');

subplot(2,1,2);
plot(f, abs(H_rc));
title('Raised Cosine脉冲的频率响应幅度');
xlabel('频率 (Hz)');
ylabel('幅度');
