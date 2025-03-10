clc;clear;close all;

% GFDM Demodulator
no_of_bin_data = 2^16;  % 二进制数据的数量设置为2的16次方
binary_source = randsrc(1,no_of_bin_data,[0 1]);  % 生成随机二进制数据源
no_of_subcarriers = 128;  % 活动子载波数量设置为128
no_of_total_subcarriers = 256;  % 总的子载波数量设置为256
len_of_CP=8;  % 循环前缀的长度
no_of_subsymbols = 64;  % 设置子符号数量
no_of_blocks = len_of_mapped_signal/(no_of_subcarriers*no_of_subsymbols);  % 计算块的数量

% 总的符号数
total_no_of_subsymbols=no_of_subsymbols;

% 脉冲成型滤波器
g = RRCPulseShape(total_no_of_subsymbols,no_of_total_subcarriers);
% 接收信号
x=receiver;
% 将接收信号重塑为矩阵形式，按块划分
y = reshape(x,[],no_of_blocks);

% 去除接收信号中的循环前缀
y = y(len_of_CP+1:end,:);

% 初始化三维数组 IQ，用于存储中间处理的信号
IQ = zeros(no_of_total_subcarriers,no_of_subsymbols,no_of_blocks);

% 初始化三维数组 D，用于存储每个子载波和子符号对应的解调结果
D = zeros(no_of_total_subcarriers,no_of_subsymbols,no_of_blocks);

% 初始化二维数组 s，用于存储每个块的解调结果
s = zeros(no_of_total_subcarriers*no_of_subsymbols,no_of_blocks);

% 遍历每个数据块进行处理
for b = 1:no_of_blocks
    % 提取当前块的接收信号
    yhat = y(:,b);

    % 获取脉冲成形滤波器的频域响应的共轭，用于匹配滤波
    G = conj(fft(g));

    % 计算每个子符号对应的频域长度
    L = length(G)/no_of_subsymbols;

    % 对当前块的接收信号进行快速傅里叶变换（FFT）
    xhat = fft(yhat);

    % 初始化二维数组 Dhat，用于存储每个子载波的解调结果
    Dhat = zeros(no_of_total_subcarriers,no_of_subsymbols);

    % 遍历每个子载波进行处理
    for k = 1:no_of_total_subcarriers
        % 对频域信号进行循环移位，分离出当前子载波的信号
        carrier = circshift(xhat,ceil(L*no_of_subsymbols/2) - no_of_subsymbols*(k-1));

        % 对分离出的子载波信号进行 FFT 移位，使其零频居中
        carrier = fftshift(carrier(1:L*no_of_subsymbols));

        % 进行匹配滤波，将子载波信号与滤波器频域响应共轭相乘
        carrierMatched = carrier .* G;

        % 对匹配滤波后的信号进行逆快速傅里叶变换（IFFT）
        % 并对每个子符号进行求和平均
        dhat = ifft(sum(reshape(carrierMatched,no_of_subsymbols,L),2)/L);

        % 将解调结果存储到 D 数组中
        D(k,:,b) = dhat;
    end

    % 将三维数组 D 中当前块的数据重塑为一维向量
    s(:,b) = reshape(D(:,:,b), numel(D(:,:,b)),1);

    % 将一维向量 s 重塑为二维矩阵 s_
    s_(:,:,b) = reshape(s(:,b),no_of_total_subcarriers,no_of_subsymbols);

    % 去除空子载波，提取有效子载波的数据
    sf(:,:,b) = s_(51:no_of_subcarriers+50,:,b);

    % 将去除空子载波后的二维矩阵 sf 重塑为一维向量
    s2(:,b) = reshape(sf(:,:,b),[],1);
end

% 将所有块的一维向量 s2 拼接成一个一维向量 s3
s3 = reshape(s2,1,[]);

% 对一维向量 s3 进行 QPSK 解调，得到二进制比特序列
sk = qamdemod(s3, 2^2,'gray','OutputType','bit','UnitAveragePower',true);

% 将解调后的二进制比特序列重塑为一维向量
sx = reshape(sk,1,[]);

% 比较解调后的二进制比特序列与原始二进制数据源
% 统计错误比特的数量
errors(si) = length(find(sx' - sen'));


% 计算不同信噪比下的误比特率（BER）
ber = errors/length(sx);




function res = RRCPulseShape(total_no_of_subsymbols,no_of_total_subcarriers)
% 该方法用于生成升余弦平方根（RRC）脉冲成形滤波器的脉冲
M = total_no_of_subsymbols;  % 获取总的子符号数量
K = no_of_total_subcarriers;  % 获取总的子载波数量
a = 0.1;  % 升余弦滤波器的滚降因子设置为0.1
% 生成时间向量
t = linspace(-M/2, M/2, M*K+1);
t = t(1:end-1);
t = t';
% 根据升余弦平方根滤波器公式计算脉冲值
g = (sin(pi*t*(1-a))+4*a.*t.*cos(pi*t*(1+a)))./(pi.*t.*(1-(4*a*t).^2));
% 处理t = 0的特殊情况
g(find(t==0)) = 1-a+4*a/pi;
% 处理|t| = 1/(4*a)的特殊情况
g(find(abs(t) == 1/(4*a))) = a/sqrt(2)*((1+2/pi)*sin(pi/(4*a))+(1-2/pi)*cos(pi/(4*a)));

g = fftshift(g);  % 将脉冲进行FFT移位
res = g / sqrt(sum(g.*g));  % 对脉冲进行归一化处理
% filter_pulse = res;  % 存储归一化后的脉冲
end