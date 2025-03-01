% 定义常量
c = 299792458; % 光速，单位：米/秒
f = 0.1e-9; % 光纤色散系数，单位：米/平方
lmbd = 1550e-9; % 光波长，单位：米
Rs = 10e9; % 符号速率，单位：符号/秒

% 设置SNR值
snr_values = [1:10]; % 假设SNR值从1到10

% 计算每个SNR值的OSNR
osnr_values = zeros(length(snr_values), 1);
for i = 1:length(snr_values)
    snr = snr_values(i);
    Bref = c * f / (lmbd.^2);
    osnr_values(i) = Rs * snr / (2 * Bref);
end

% 输出计算结果
disp('计算得到的OSNR值：');
disp(osnr_values);

% % 假设的理论OSNR值（这里仅用于验证，实际计算中不需要）
% theoretical_osnr_values = [1:10]; % 假设的理论OSNR值从1到10
% 
% % 比较计算结果与理论值
% error_magnitude = abs(osnr_values - theoretical_osnr_values);
% disp('计算误差：');
% disp(error_magnitude);
