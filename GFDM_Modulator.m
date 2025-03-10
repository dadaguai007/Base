clc;clear;close all;

% GFDM 调制

% % 基本参数
% no_of_subcarriers;  % 活动子载波的数量
% total_no_of_subsymbols;  % 总的子符号数量
% no_of_subsymbols;  % 子符号的数量
% no_of_blocks;  % 块的数量
% no_of_total_subcarriers;  % 总的子载波数量
% no_of_bin_data;  % 二进制数据的数量
% len_of_mapped_signal;  % 映射后信号的长度
% binary_source;  % 二进制数据源
% maped_signal;  % 映射后的信号
% no_of_bits_per_symbol;  % 每个符号的比特数



no_of_bin_data = 2^16;  % 二进制数据的数量设置为2的16次方
binary_source = randsrc(1,no_of_bin_data,[0 1]);  % 生成随机二进制数据源
no_of_subcarriers = 128;  % 活动子载波数量设置为128
no_of_total_subcarriers = 256;  % 总的子载波数量设置为256



no_of_bits_per_symbol = 4;  % 16QAM调制每个符号携带4比特
% 使用qammod函数进行16QAM调制
maped_signal_ = qammod(binary_source', 16, 'InputType','bit','UnitAveragePower',true);
maped_signal = maped_signal_;  % 存储映射后的信号

len_of_mapped_signal = length(maped_signal_);  % 记录映射后信号的长度




cp_len=8;  % 循环前缀的长度
maped_sig = maped_signal;  % 获取映射后的信号

% GFDM Modulator

no_of_subsymbols = 64;  % 设置子符号数量
no_of_blocks = len_of_mapped_signal/(no_of_subcarriers*no_of_subsymbols);  % 计算块的数量

block_len = no_of_subsymbols * no_of_subcarriers;  % 计算块的长度
total_no_of_subsymbols = no_of_subsymbols;  % 记录总的子符号数量（与输入的子符号数量相同）
total_block_len = no_of_total_subcarriers * total_no_of_subsymbols;  % 计算总的块长度
% 初始化GFDM信号向量
GFDM_signal = zeros((no_of_blocks*total_block_len)+cp_len,1);

for i = 0:no_of_blocks-1
    % 提取当前块的信号
    sig = maped_sig(block_len*i+1:end-(block_len*(no_of_blocks - 1 - i)));
    % 对信号进行重塑
    sig_reshaped = reshape(sig, no_of_subcarriers, no_of_subsymbols);
    % 初始化总的信号矩阵
    sig_total = zeros(no_of_total_subcarriers,total_no_of_subsymbols);
    % 将重塑后的信号放置在指定位置
    sig_total(51:no_of_subcarriers+50,:) = sig_reshaped;
    % 生成升余弦平方根脉冲成形滤波器的脉冲
    filter_value = RRCPulseShape(total_no_of_subsymbols,no_of_total_subcarriers);
    % 对信号进行逆快速傅里叶变换并重复
    sig_rep = repmat(no_of_total_subcarriers*ifft(sig_total),total_no_of_subsymbols,1);

    % 初始化  求和信号向量
    sig_sum = zeros(total_block_len, 1);

    for j = 1:no_of_subsymbols
        % 对每个子符号进行脉冲成形
        sym = sig_rep(:,j).*filter_value;
        % 对信号进行循环移位
        sym = circshift(sym, no_of_total_subcarriers*(j-1));
        % 累加信号
        sig_sum = sig_sum + sym;
    end
    % 添加循环前缀
    sig_cp = [sig_sum(end-cp_len+1:end); sig_sum];
    % 将当前块的信号添加到总的GFDM信号中
    GFDM_signal(i*(total_block_len+cp_len) + (1:total_block_len + cp_len)) = sig_cp;
end
% 将GFDM信号转换为一维向量
transmitted_signal = GFDM_signal';
% 对发射信号进行归一化处理
transmitted_signal = transmitted_signal/std(transmitted_signal);





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