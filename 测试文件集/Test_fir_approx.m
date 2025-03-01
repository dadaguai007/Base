% Test_fir_approx
clc;clear;
% 定义一个简单的通道频率响应，例如一个低通滤波器
Hch = @(f) double(f < 0.2); % 一个理想的低通滤波器，截止频率为0.2

% 设置采样频率
fs = 1; % 采样频率为1

% 设置能量分数
Efraction = 0.95; % 我们希望FIR滤波器包含至少95%的能量

% 调用fir_approx函数
[Ntaps, hfir] = fir_approx(Hch, fs, Efraction);

% 显示结果
disp(['FIR滤波器的抽头数: ', num2str(Ntaps)]);
disp('FIR滤波器的系数:');
disp(hfir);
