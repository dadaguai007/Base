%FBMC Demodulator

% The example implements a basic FBMC demodulator and measures the BER for
% the chosen configuration. The processing
% includes matched filtering followed by OQAM separation to form the
% received data symbols. These are demapped to bits and the resultant bit
% error rate is determined.

% 定义FFT点数，用于频域处理
N = 1024;           % Number of FFT points
% 定义双边保护带的数量，避免频谱泄漏
numGuards = 212;         % Guard bands on both sides
% 重叠符号数量，影响时域符号的重叠特性
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
% MQAM调制阶数
M = 4;
% 每个符号携带的比特数，根据调制阶数计算
k = log2(M);   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM

% 直接指定OFDM符号数量
numOFDMSymbols = 10;

%   Initialize arrays
% 计算携带数据的子载波数量
L = N - 2*numGuards;  % Number of subcarriers with data = QAM symbols per OFDM symbol
% K倍的FFT点数
KN = K*N;
% K倍的携带数据的子载波数量
KL = K*L;

% We load the results of our experiment as well as the input data to
% compute the BER
% 加载发射信号数据
load('SignalTransmitted.mat')
% 加载输入比特数据，用于计算误比特率
load('InputData.mat')


% QAM demodulator
% 创建矩形QAM解调器对象，用于将QAM符号解调为比特序列
qamDemod = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^k, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

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
        % 如果K不是2、3、4，结束程序
        return
end

% Build symmetric filter
% 构建对称的原型滤波器系数
Hk = [fliplr(HkOneSided) 1 HkOneSided];


% 接收信号
SignalRecibida=SignalRec;

% Convert from real to complex
% 将一维的接收信号重新排列为二维矩阵，准备转换为复信号
fbmc_real = reshape(SignalRecibida, K*N, 2*numOFDMSymbols);
% 初始化复信号矩阵
txSigAll = zeros(K*N, numOFDMSymbols);
% 遍历每个OFDM符号，将实信号转换为复信号
for symIdx = 1:numOFDMSymbols
    txSigAll(:, symIdx) = fbmc_real(:, (2*symIdx) - 1) + 1i*fbmc_real(:, 2*symIdx);
end

% Process symbol-wise
% 逐符号处理接收信号
for symIdx = 1:numOFDMSymbols
    % 提取当前符号的接收信号
    rxSig = txSigAll(:, symIdx);

    % Perform FFT
    % 对接收信号进行FFT变换，将信号从时域转换到频域
    rxf = fft(fftshift(rxSig));

    % Matched filtering with prototype filter
    % 使用原型滤波器对频域信号进行匹配滤波
    rxfmf = filter(Hk, 1, rxf);
    % Remove K - 1 delay elements
    % 移除滤波器引入的K - 1个延迟元素
    rxfmf = [rxfmf(K:end); zeros(K - 1, 1)];
    % Remove guards
    % 移除保护带
    rxfmfg = rxfmf(numGuards*K + 1:end - numGuards*K);

    % OQAM post-processing
    %  Downsample by 2K, extract real and imaginary parts
    % OQAM后处理，按2K进行下采样，提取实部和虚部
    if rem(symIdx, 2)
        % 奇数符号：虚部在实部之后K个样本
        r1 = real(rxfmfg(1:2*K:end));
        r2 = imag(rxfmfg(K + 1:2*K:end));
        rcomb = complex(r1, r2);
    else
        % 偶数符号：实部在虚部之后K个样本
        r1 = imag(rxfmfg(1:2*K:end));
        r2 = real(rxfmfg(K + 1:2*K:end));
        rcomb = complex(r2, r1);
    end
    %  Normalize by the upsampling factor
    % 按上采样因子K进行归一化
    rcomb = (1/K)*rcomb;

    % Demapper: Perform hard decision
    % 使用QAM解调器进行硬判决，将QAM符号解调为比特序列
    rxBits(:, symIdx) = qamDemod(rcomb);
end

% Measure BER with appropriate delay
% 创建误比特率计算对象
BER = comm.ErrorRate;
% 设置接收延迟，确保输入数据和接收比特序列对齐
BER.ReceiveDelay = k*KL;
% 计算误比特率
ber = BER(inpData(:), rxBits(:));

% Display Bit error
% 显示当前衰减和K值下的误比特率
disp(['FBMC Reception for K = ' num2str(K) ', BER = ' num2str(ber(1)) ])

% 保存误比特率矩阵到文件
save(nameMatrixBERatt, 'matrixBER')