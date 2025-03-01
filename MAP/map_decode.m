% MAP
function [L_ext, decoded_bits]=map_decode(sys_symbols, non_sys_symbols, trellis, noise_variance, L_int, mode)
% Parameters
% ----------
% sys_symbols : 1D ndarray
% Received symbols corresponding to
% the systematic (first output) bits in
% the codeword.
%
% non_sys_symbols : 1D ndarray
% Received symbols corresponding to the non-systematic
% (second output) bits in the codeword.
%
% trellis : Trellis object
% Trellis representation of the convolutional code.
%
% noise_variance : float
% Variance (power) of the AWGN channel.
%
% L_int : 1D ndarray
% Array representing the initial intrinsic
% information for all received
% symbols.
%
% Typically all zeros,
% corresponding to equal prior
% probabilities of bits 0 and 1.
%
% mode : str{'decode', 'compute'}, optional
% The mode in which the MAP decoder is used.
% 'decode' mode returns the decoded bits
%
% along with the extrinsic information.
% 'compute' mode returns only the
% extrinsic information.

if nargin<6
    mode='decode';
end
% 状态转移表参数
k = trellis.k;
n = trellis.n;
number_states = trellis.number_states;
number_inputs = trellis.number_inputs;

% 输入长度
msg_length = length(sys_symbols);

%  Initialize forward state metrics (alpha),初始化前馈状态
f_state_metrics = zeros(number_states, 2);
f_state_metrics(1,1) = 1;


% # Initialize backward state metrics (beta),初始化回溯状态
b_state_metrics = zeros(number_states, msg_length+1);
b_state_metrics(:,msg_length) = 1;

%  Initialize branch transition probabilities (gamma) 度量矩阵初始化
branch_probs = zeros(number_inputs, number_states, msg_length+1);

% 后验概率；
app = zeros(number_inputs,1);
% LLR初始化
% lappr = 0; % LLR

% decoded_bits = zeros(msg_length, 'int')
L_ext = zeros(msg_length,1);
% 先验概率
priors = zeros(2, msg_length);
% 先验概率计算
priors(1,:) = 1/(1 + exp(L_int)); % 概率密度公式
priors(2,:) = 1 - priors(1,:);  % 两者相加为1

% # Backward recursion 回溯矩阵
b_state_metrics=backward_recursion(trellis, msg_length, noise_variance, sys_symbols,non_sys_symbols, branch_probs, priors, b_state_metrics);
%
% # Forward recursion  前馈递进
[L_ext,decoded_bits]=forward_recursion_decoding(trellis, mode, msg_length, b_state_metrics,f_state_metrics, branch_probs, app, L_int,priors, L_ext, decoded_bits);


end