clc;clear;close all;

% 说明：
% FBMC 会对每个经过调制的子载波信号进行滤波处理，目的是降低带外辐射。
% 在多载波通信系统中，不同子载波在频谱上是相邻的，如果子载波信号的带外辐射过大，
% 就会对相邻子载波产生干扰，影响整个系统的性能。
% 通过滤波可以使子载波信号的能量更集中在其分配的频带内，减少对其他子载波的干扰。

% 滤波器的特性由重叠因子 K 来表征。
% K 表示在时域中相互重叠的多载波符号的数量。
% 在 FBMC 系统中，符号在时域上是重叠的，这种重叠特性有助于提高频谱效率和抗干扰能力。
% 不同的 K 值会影响滤波器的性能和系统的整体表现。

% 原型滤波器的阶数可以选择为 2*K - 1，其中 K 的取值可以是 2、3 或 4。

%  FBMC 实现采用了频率扩展技术。使用长度为 N*K 的快速傅里叶逆变换（IFFT），
% 并且符号之间存在 N/2 的延迟，其中 N 是子载波的数量。
% 频率扩展可以将信号的能量在频域上更均匀地分布，提高系统的抗干扰能力。
% 而符号之间的延迟设计是为了实现时域上的符号重叠，使得相邻符号之间可以相互协作，进一步提高系统性能。
% 这种设计选择使得对 FBMC 系统的分析变得容易，并且便于与其他调制方法进行比较。

% 采用了固定的 K 值、特定的原型滤波器阶数、频率扩展和符号延迟等设计。


% FBMC基础之一是使用OOQAM调制

%FBMC Modulator
% 设置FFT点数，即频域处理的点数
N = 1024;           % Number of FFT points
% 双边保护带的数量，用于避免频谱泄漏和邻道干扰
numGuards = 212;         % Guard bands on both sides
% 重叠符号的数量，FBMC中时域符号的重叠因子，取值可以是2、3或4
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
% MQAM调制阶数
M = 4; 
% 每个符号携带的比特数，根据调制阶数计算得出
k = log2(M);   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM

% 直接设定OFDM符号的数量
numOFDMSymbols = 10;



% Prototype filter
% 根据重叠因子K选择原型滤波器的单边系数
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        % 如果K不是2、3或4，直接返回，结束程序
        return
end

% Build symmetric filter
% 构建对称的原型滤波器系数，将单边系数翻转后拼接
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% Transmit-end processing
%   Initialize arrays
% 计算携带数据的子载波数量，即每个OFDM符号中的QAM符号数量
L = N - 2*numGuards;  % Number of subcarriers with data = QAM symbols per OFDM symbol
% K倍的FFT点数
KN = K*N;
% K倍的携带数据的子载波数量
KL = K*L;
% 初始化携带数据的子载波数组
dataSubCar = zeros(L, 1);
% 初始化上采样后的子载波数据数组
dataSubCarUp = zeros(KL, 1);

% 每个符号的比特数，除以2是因为采用了偏移正交幅度调制（OQAM）
numBits = k*L/2;    % The 2 is as we will use double space for OQAM

% 初始化输入数据矩阵，存储每个符号的输入比特
inpData = zeros(numBits, numOFDMSymbols); 
% 初始化接收比特矩阵，这里暂时未在调制器中使用
rxBits = zeros(numBits, numOFDMSymbols);
% 初始化发射信号矩阵，存储每个符号的发射信号
txSigAll = complex(zeros(KN, numOFDMSymbols));
% 初始化符号缓冲区，用于存储延迟的符号
symBuf = complex(zeros(2*KN, 1));

% We configure the QAM symbol mapper
% 创建矩形QAM调制器对象，将输入比特映射为QAM符号
qamMapper = comm.RectangularQAMModulator(...
    'ModulationOrder', M, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% For each column we generate the data (inputData) and perform the Filter Bank
% 遍历每个OFDM符号
for symIdx = 1:numOFDMSymbols
    % Generate data and map qam
    % 生成随机比特序列作为输入数据
    inpData(:, symIdx) = randi([0 1], numBits, 1);
    % 使用QAM调制器将输入比特映射为QAM符号
    modData = qamMapper(inpData(:, symIdx));
    
    % OQAM Modulator: alternate real and imaginary parts
    % 偏移正交幅度调制（OQAM），奇数符号和偶数符号交替处理实部和虚部
    if rem(symIdx, 2) == 1     % Odd symbols
        % 奇数符号：将调制符号的实部放入奇数索引的子载波，虚部乘以j放入偶数索引的子载波
        dataSubCar(1:2:L) = real(modData);
        dataSubCar(2:2:L) = 1i*imag(modData);
    else                    % Even symbols
        % 偶数符号：将调制符号的虚部乘以j放入奇数索引的子载波，实部放入偶数索引的子载波
        dataSubCar(1:2:L) = 1i*imag(modData);
        dataSubCar(2:2:L) = real(modData);
    end

    % Upsample by K, pad with guards, and filter with the prototype filter
    % 对携带数据的子载波数据进行K倍上采样
    dataSubCarUp(1:K:end) = dataSubCar;
    % 在数据两端添加保护带
    dataBitsUpPad = [zeros(numGuards*K, 1); dataSubCarUp; zeros(numGuards*K, 1)];
    
    % 使用原型滤波器对数据进行滤波
    X1 = filter(Hk, 1, dataBitsUpPad);
    
    % Remove 1/2 filter length delay
    % 去除滤波器引入的1/2滤波器长度的延迟
    X = [X1(K:end); zeros(K - 1, 1)];

    % Compute IFFT of length KN for the transmitted symbol
    % 对处理后的数据进行长度为KN的IFFT变换，并将零频移到中心
    txSymb = fftshift(ifft(X));

    % Transmitted signal is a sum of the delayed real, imag symbols
    % 更新符号缓冲区，将前面的部分移除，后面添加零
    symBuf = [symBuf(N/2 + 1:end); complex(zeros(N/2, 1))];
    % 将当前符号叠加到符号缓冲区的相应位置
    symBuf(KN + (1:KN)) = symBuf(KN + (1:KN)) + txSymb;

    % Compute power spectral density (PSD) for the plots later on
    % 提取当前符号的信号
    currSym = complex(symBuf(1:KN));
    
    % Store transmitted signals for all symbols
    % 将当前符号的发射信号存储到发射信号矩阵中
    txSigAll(:, symIdx) = currSym;
end

%将生成的信号转换为适合实时处理的形式。
% 在实际通信系统中，实时处理对信号的格式和特性有特定要求，所以需要进行这样的转换

% 将传统的FBMC复数信号转换为实信号，便于实时处理
fbmc_real = zeros(K*N, 2*numOFDMSymbols);
for symIdx = 1:numOFDMSymbols
    % 偶数列存储虚部
    fbmc_real(:,(symIdx*2)) = imag(txSigAll(:,symIdx));
    % 奇数列存储实部
    fbmc_real(:,(symIdx*2)-1) = real(txSigAll(:,symIdx));
end

% Plot original power spectral density
% 初始化FBMC信号和OFDM信号的功率谱密度总和数组
sumFBMCSpec = zeros(KN*2, 1);
sumOFDMSpec = zeros(N*2, 1);
% 遍历每个OFDM符号
for symIdx = 1:numOFDMSymbols
    % 提取当前符号的实信号
    currentSymbol = fbmc_real(:,symIdx*2);
    % 使用周期图法计算当前符号的功率谱密度
    [specFBMC, fFBMC] = periodogram(complex(currentSymbol), hann(KN, 'periodic'), KN*2, 1);
    % 累加功率谱密度
    sumFBMCSpec = sumFBMCSpec + specFBMC;
end
% 对FBMC信号的功率谱密度进行归一化处理
sumFBMCSpec = sumFBMCSpec/mean(sumFBMCSpec(1 + K + 2*numGuards*K:end - 2*numGuards*K - K));
% 绘制FBMC信号的功率谱密度曲线
plot(fFBMC - 0.5, 10*log10(sumFBMCSpec));
% 显示网格线
grid on
% 设置坐标轴范围
axis([-0.5 0.5 -180 10]);
% 设置x轴标签
xlabel('Normalized frequency'); 
% 设置y轴标签
ylabel('PSD (dBW/Hz)')
% 设置图的标题
title(['Real FBMC, K = ' num2str(K) ' overlapped symbols'])
% 设置图形窗口的位置
set(gcf, 'Position', figposition([15 50 30 30]));

%Paralel to serie and clocks
% 将实信号矩阵转换为一维向量，实现并行转串行
SignalTransmitted = reshape(fbmc_real, K*N*2*numOFDMSymbols, 1);   


% 保存发射信号到文件
save('SignalTransmitted.mat', 'SignalTransmitted')
% 保存输入数据到文件
save('InputData.mat', 'inpData')