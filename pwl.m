function [y,e,w] = pwl(sig,train,ntaps,step,thVec)
%它接收五个输入参数：sig（输入信号）、train（训练信号或目标信号）、
% ntaps（滤波器的抽头数）、step（步长，用于调整滤波器权重）、
% thVec（阈值向量，用于确定分解信号的阈值）

% signal decomposition
N = length(thVec)+1;
%将输入信号sig分解为不同的波形分量，将用于初始化滤波器的权重
M = cal_pwl_decomposition_table(sig,thVec);
% filter initiation
% 初始化权重，中间抽头初始化为1
w = zeros(ntaps,N);
w(ceil(ntaps/2),:) = 1;
% start
for idx = 1:length(sig)-ntaps+1
    % 直接按照矩阵进行选取信号
    x_window = M(idx:idx+ntaps-1,:);
    %均衡信号并进行加权
    y(idx) = sum(sum(w.*x_window));
    %更新系数
    if idx+ceil(ntaps/2)-1 < length(train)
        e(idx) = train(idx+ceil(ntaps/2)-1) - y(idx);
    else
        e(idx) = sign(y(idx)) - y(idx);
    end
    w = w+step*e(idx)*x_window;
end

% 如何调用：
% [y,e] = pwl_bo(noisySig,ref,9,0.0003,[-0.7661,0,0.7661]);
% ref = modData./(M-1);