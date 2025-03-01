% Test_中值滤波
% 创建一个简单的输入信号，例如一个带有噪声的正弦波
t = 0:0.01:1; % 时间向量，从0到1秒，步长为0.01秒
f = 5; % 正弦波的频率为5Hz
x = sin(2*pi*f*t); % 创建纯净的正弦波
x_noisy = x + 0.5*randn(size(x)); % 添加高斯白噪声

% 设置滤波次数和滤波器阶数
m = 3; % 滤波3次
order = 5; % 使用5阶中值滤波器

% 调用multimidfilter函数
y = multimidfilter(x_noisy, m, order);

% 绘制原始信号、噪声信号和多重中值滤波后的信号
figure;
subplot(3,1,1);
plot(t, x);
title('原始信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(3,1,2);
plot(t, x_noisy);
title('噪声信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(3,1,3);
plot(t, y);
title('多重中值滤波后的信号');
xlabel('时间 (s)');
ylabel('幅度');
