% Test_H_imdd
% 设置参数
f = 0:10:1e6; % 频率向量，从0到10000 Hz
L = 40e3; % 光纤长度，40 km
D = 17; % 色散参数，假设为10 ps/nm/km
alpha = 0.2; % 色散参数，假设为0.5
type = 'small signal'; % 类型为 'large signal'

% 调用函数
[Hf] = Himdd(f, L, D, alpha, type);

% 显示结果
disp('Frequency response:');
disp(Hf);
figure;
plot(Hf)
