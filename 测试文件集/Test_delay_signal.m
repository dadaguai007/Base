% Test_delay_signal

% 创建一个测试信号：一个频率为 5 Hz 的正弦波，采样率为 1000 Hz，持续 1 秒
fs = 1000; % 采样率
t = 0:1/fs:1-1/fs; % 时间向量
f = 5; % 正弦波频率
x = sin(2*pi*f*t); % 创建正弦波信号


% 绘制原始信号
figure;
plot(t, x);
title('原始信号');
xlabel('时间 (s)');
ylabel('幅度');

% 应用整数延迟
delay_int = 25; % 延迟 50 个样本
y_int = delay_signal(x, delay_int);

% 绘制整数延迟后的信号
figure;
plot(t, y_int);
title('整数延迟后的信号');
xlabel('时间 (s)');
ylabel('幅度');

% 应用非整数延迟
delay_float = 50.5; % 延迟 50.5 个样本
y_float = delay_signal(x, delay_float);

% 绘制非整数延迟后的信号
figure;
plot(t, y_float);
title('非整数延迟后的信号');
xlabel('时间 (s)');
ylabel('幅度');

% 比较原始信号和延迟后的信号
figure;
plot(t, [x;y_int;y_float], 'linewidth', 2);
legend('原始信号', '整数延迟', '非整数延迟');
title('原始信号与延迟后的信号对比');
xlabel('时间 (s)');
ylabel('幅度');
