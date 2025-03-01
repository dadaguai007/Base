clc;clear;
% chirp linear simulation

Fa = 200;  % 采样频率
Ta = 1 / Fa;  % 采样周期
B = 100;

d = 1600;
t = 0:Ta/16:0.5;  % 生成时间向量，采样时间从0到0.5秒，间隔为Ta/16

f0 = -100;
f1 = 100;
t1 = max(t);
x = chirp(t, f0, t1, f1, 'linear');  % 生成线性啁啾信号

% 下采样
xa = x(1:16:end);
ta = t(1:16:end);

% 绘制波形
figure;
plot(t, x, '-', 'DisplayName', 'y(t)');
hold on;
plot(ta, xa, 'ko', 'MarkerSize', 4, 'DisplayName', 'x[kTa]');
xlabel('时间 (s)');
grid on;
legend;
xlim([min(t), max(t)]);


%
% T=10e-5;                  % 脉冲宽度100微秒
% B=1e6;                    % 带宽1MHz
% k=B/T;                    % 调频斜率
% fs=5e6;
% N=fs*T;
% t=linspace(-T/2,T/2,N); 
% s=exp(1i*k*pi*t.^2);      % 线性调频（LFM）信号
% figure;
% plot(t,s);                  
% title('线性调频(LFM)信号');
% xlabel('t/s');ylabel('幅度');%单个脉冲