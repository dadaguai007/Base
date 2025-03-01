function   outputBits=conv_encoder(inputBits,constraintLength)
%conv_encoder卷积编码器
% 参数设置[2,1,3]编码器
% inputBits = [1 0 1 1 0 1]; % 输入比特序列
% constraintLength = 3; % 约束长度
numInputBits = length(inputBits); % 输入比特序列的长度
numOutputBits = numInputBits * 2; % 输出比特序列的长度

% 初始化编码器状态
state = zeros(1, constraintLength - 1); % 初始状态全为0

% 初始化编码器输出
outputBits = zeros(1, numOutputBits);

% 编码操作
%状态也是由V2V1表示的
for i = 1:numInputBits
    % 计算编码输出
    %输入bit和状态2的异或，对应输出v1
    output1 = mod(inputBits(i) + state(2), 2);
    %输入bit和状态1，状态2异或，对应输出v0
    output2 = mod(inputBits(i) + state(1) + state(2), 2);

    % 更新状态
    state = [inputBits(i), state(1), state(2)];

    % 存储输出
    %输出的格式v1v0，
    outputBits(2*i - 1) = output1;
    outputBits(2*i) = output2;
end

% disp('输入比特序列：');
% disp(inputBits);
% disp('编码器输出比特序列：');
% disp(outputBits);
end