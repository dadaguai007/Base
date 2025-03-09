clc;clear;close all;
% 恢复信号为复数形式，并转换为矩阵形式；进行相应的均衡；转为一维向量，进行解码

% UFMC Demodulator
% 设置FFT点数，用于后续的频域处理，要与调制端保持一致
numFFT = 1024;        % number of FFT points
% 设置每个子带的大小，该值必须大于1，与调制端一致
subbandSize = 30;    % must be > 1
% 设置子带的数量，需满足子带数量乘以子带大小小于等于FFT点数，与调制端一致
numSubbands = 20;    % numSubbands*subbandSize <= numFFT
% 设置子带的偏移量，用于将信号放置在频谱的中心位置，与调制端一致
subbandOffset = 212; % numFFT/2-subbandSize*numSubbands/2 for band center
% 设置UFMC符号的数量，决定了要处理的符号个数，与调制端一致
numUFMCSymbols = 10;
% 加载保存的发射信号数据
load('SignalTransmitted.mat')

% Dolph - Chebyshev窗设计参数
% 滤波器的长度，类似于循环前缀的长度，与调制端一致
filterLen = 43;      % similar to cyclic prefix length
% 旁瓣衰减，单位为dB，用于控制滤波器的旁瓣特性，与调制端一致
slobeAtten = 40;     % sidelobe attenuation, dB
% MQAM调制阶数，这里设置为4，表示采用4QAM调制，与调制端一致
M = 4;               % MQAM: 4QAM, 16QAM, 64QAM, 256QAM
% 每个符号携带的比特数，通过对数运算得到，与调制端一致
k = log2(M);

% Design window with specified attenuation
% 使用chebwin函数设计Dolph - Chebyshev窗作为原型滤波器，与调制端一致
prototypeFilter = chebwin(filterLen, slobeAtten);



% 接收信号（待使用者进行更改）
SignalRecibida=SignalRec;


% Serie to paralel and real to complex
% 将一维的接收信号重新排列成二维矩阵，实现并行化
UFMC_real = reshape(SignalRecibida, numFFT+filterLen-1, 2*numUFMCSymbols);
% 初始化一个矩阵，用于存储复数形式的接收信号
txSigAll = zeros(numFFT+filterLen-1,numUFMCSymbols);
% 遍历每个UFMC符号
for symIdx = 1:numUFMCSymbols
    % 将实部和虚部组合成复数信号
    txSigAll(:,(symIdx)) = UFMC_real(:,(2*symIdx)-1)+1i*UFMC_real(:,(2*symIdx));
end
% 初始化一个矩阵，用于存储解调后的比特序列
rxBits = zeros(k*subbandSize*numSubbands,numUFMCSymbols);

% 遍历每个UFMC符号进行解调处理
for symIdx = 1:numUFMCSymbols
    % 提取当前UFMC符号的接收信号
    yRx = txSigAll(:,(symIdx));

    % 将接收信号向量填充到FFT长度的两倍，这里没有采用加窗或额外滤波操作
    yRxPadded = [yRx; zeros(2*numFFT - numel(yRx),1)];

    % Perform FFT and downsample by 2
    % 对填充后的信号进行FFT变换并进行降采样，采样间隔为2
    RxSymbols2x = fftshift(fft(yRxPadded));
    RxSymbols = RxSymbols2x(1:2:end);

    % Select data subcarriers
    % 从FFT结果中选择数据子载波对应的符号
    dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));


    % Use zero - forcing equalizer after OFDM demodulation
    % 在OFDM解调后使用迫零均衡器
    rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen - 1)'/numFFT); ...
        zeros(numFFT - filterLen,1)];
    % 对滤波器进行FFT变换并移到中心位置
    prototypeFilterFreq = fftshift(fft(rxf));
    % 计算原型滤波器在数据子载波上的逆频率响应（选取每个载波上对应的滤波器响应）
    prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2 - subbandSize/2+(1:subbandSize));

    % Equalize per subband - undo the filter distortion
    % 将数据子载波符号重新排列成矩阵形式
    dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
    % 对每个子带进行均衡，消除滤波器引入的失真
    EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
    % 将均衡后的矩阵转换为一维向量
    EqualizedRxSymbols = EqualizedRxSymbolsMat(:);


    % Demapping and BER computation
    % 创建一个矩形QAM解调器对象，用于将QAM符号解调为比特序列
    qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
        2^k, 'BitOutput', true, ...
        'NormalizationMethod', 'Average power');
    % 创建一个误比特率计算对象
    BER = comm.ErrorRate;

    % Perform hard decision and measure errors
    % 对均衡后的符号进行硬判决，得到解调后的比特序列
    rxBits(:,(symIdx)) = qamDemod(EqualizedRxSymbols);
end

% 再次加载输入数据
load('inputdata.mat')
% 计算误比特率
ber = BER(inputData(:), rxBits(:));

% 显示当前接收情况下的误比特率
disp(['UFMC Reception, BER = ' num2str(ber(1))]);

