clc;clear
% 生成随机电信号
fb=40e3;
fb2=60e3;
fb3=80e3;
fb4=100e3;
fb5=20e3;
Fs = 2e9;  % 采样率
% Fs=250e6; %机器自动识别采样率
N=100000;
t = 0:(1/Fs):N*(1/Fs)-1/Fs;  % 时间向量
f= Fs * (-0.5:1/N:0.5-1/N);
x=cos(2*pi*fb*t)+cos(2*pi*fb2*t)+cos(2*pi*fb3*t)+cos(2*pi*fb4*t)+cos(2*pi*fb5*t);
y=sin(2*pi*fb*t)+sin(2*pi*fb2*t)+sin(2*pi*fb3*t)+sin(2*pi*fb4*t)+sin(2*pi*fb5*t);
figure;hold on;
plot(x)
plot(y)
T=N/250e6;
disp(T)

% 将时间和信号合并
signal_x = [t', x'];
signal_y = [t', y'];
% 使用csv进行存储
% csvwrite('signal_x.csv',signal_x);
% csvwrite('signal_y.csv',signal_y);