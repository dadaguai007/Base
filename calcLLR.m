function LLRs = calcLLR(rxSymb, sigma2, constSymb, bitMap, px)
    % LLR calculation 

    M = length(constSymb);
    b = log2(M);
    
    LLRs = zeros(1, length(rxSymb) * b);

    for i = 1:length(rxSymb)
        prob = exp(-abs(rxSymb(i) - constSymb).^2 / sigma2) .* px;
        
        for indBit = 1:b
            p0 = sum(prob(bitMap(:, indBit) == 0));
            p1 = sum(prob(bitMap(:, indBit) == 1));
            % 索引问题,应该从第0索引开始出发
            LLRs((i-1)* b + indBit) = log(p0) - log(p1);
        end
    end
end
%这是一个LLR的计算公式，表示第i个接收到的符号对应的第indBit个比特的LLR值。其中，p0是满足bitMap[:, indBit]等于0的符号的概率和，p1是满足bitMap[:, indBit]等于1的符号的概率和。LLR值是p0的对数与p1的对数之差。

%该函数用于计算经过循环均匀加性白高斯噪声信道后的符号序列的LLR。LLR是信道输出对符号概率的对数比，用于在信号解调过程中进行符号的软判决。函数接受五个参数，分别是：

% rxSymb：接收到的符号序列，类型为numpy数组；
% σ2：噪声方差，类型为标量；
% constSymb：星座符号，类型为形状为(M, 1)的numpy数组；
% px：先验符号概率，类型为形状为(M, 1)的numpy数组；
% bitMap：比特到符号的映射，类型为形状为(M, log2(M))的numpy数组。
% 函数返回计算得到的LLR序列。
% 在函数内部，首先获取星座符号的个数M，并计算出比特数b=log2(M)。
% 然后创建一个形状为(len(rxSymb) * b)的全零数组LLRs，用于存储计算得到的LLR值。
% 接下来，使用双重循环对每个符号进行计算。外层循环遍历接收到的符号序列rxSymb，内层循环遍历比特数b。
% 对于每个符号，首先计算其到星座符号间的幅度的负绝对值的指数，然后乘以先验符号概率px。
% 然后根据比特到符号的映射bitMap，将符合条件的项进行求和，得到p0和p1。
% 最后，将LLR值计算为对p0取对数减去对p1取对数，并将结果存储在LLRs数组的相应位置。
% 最后，返回计算得到的LLR序列LLRs。