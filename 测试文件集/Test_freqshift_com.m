% Test_freqshift
clc;clear;
% 定义一个时间向量
t = 0:0.001:1; % 从0到1秒，步长为0.001秒
Fs=1/0.001;
% 定义一个简单的输入信号，例如一个正弦波
f0 = 5; % 原始频率为5Hz
x = sin(2*pi*f0*t);

% 设置频率偏移量
fshift = 7; % 频率偏移量为2Hz

% 调用freqshift函数
xshift = freqshift(x, t, fshift);
% 调用comp_freq_offset函数
comp_signal = comp_freq_offset(xshift, fshift, Fs);
% 绘制原始信号和频率移位后的信号
figure;
subplot(3,1,1);
plot(t, x);
title('原始信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(3,1,2);
plot(t, xshift);
title('频率移位后的信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(3,1,3);
plot(t, comp_signal)

