clc;clear;close all;
% UFMC Modulator

% UFMC , 设置子带信号，每个子带信号中有N个数量，小于FFT点数
% 每个子带进行滤波操作；
% 符号数自行确定；偏移量

% 设置FFT点数，即快速傅里叶变换的点数，用于后续的频域处理
numFFT = 1024;        % number of FFT points
% 设置每个子带的大小，该值必须大于1
subbandSize = 30;     % must be > 1 
% 设置子带的数量，需满足子带数量乘以子带大小小于等于FFT点数
numSubbands = 20;     % numSubbands*subbandSize <= numFFT
% 设置子带的偏移量，通常用于将信号放置在频谱的中心位置
subbandOffset = 212;  % numFFT/2-subbandSize*numSubbands/2 for band center
% 设置UFMC符号的数量，决定了要处理的符号个数
numUFMCSymbols = 10;

% Dolph - Chebyshev窗设计参数
% 滤波器的长度，类似于循环前缀的长度
filterLen = 43;       % similar to cyclic prefix length
% 旁瓣衰减，单位为dB，用于控制滤波器的旁瓣特性
slobeAtten = 40;      % sidelobe attenuation, dB
% MQAM调制阶数，这里设置为4，表示采用4QAM调制
M = 4;                % MQAM: 4QAM, 16QAM, 64QAM, 256QAM
% 每个符号携带的比特数，通过对数运算得到
k = log2(M);   

% Design window with specified attenuation
% 使用chebwin函数设计Dolph - Chebyshev窗作为原型滤波器
prototypeFilter = chebwin(filterLen, slobeAtten);

% QAM Symbol mapper
% 创建一个矩形QAM调制器对象，用于将输入的比特序列映射为QAM符号
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^k, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit - end processing
%  Initialize arrays
% 初始化一个矩阵，用于存储每个子带的输入比特序列
inpData = zeros(k*subbandSize, numSubbands);
% 初始化一个矩阵，用于存储所有UFMC符号的输入比特序列
inputData = zeros(k*subbandSize*numSubbands,numUFMCSymbols);
% 初始化一个复数矩阵，用于存储每个UFMC符号的复数信号
UFMCcomplex = complex(zeros(numFFT+filterLen-1, numUFMCSymbols));

% 遍历每个UFMC符号
for symIdx = 1:numUFMCSymbols
    % 初始化一个复数向量，用于存储当前UFMC符号的发射信号
    txSig = complex(zeros(numFFT+filterLen-1, 1));
    % 遍历每个子带
    for bandIdx = 1:numSubbands
        % 生成随机比特序列，长度为每个子带的比特数
        bitsIn = randi([0 1], k*subbandSize, 1);
        % 使用QAM调制器将比特序列映射为QAM符号
        symbolsIn = qamMapper(bitsIn);
        % 记录当前子带的输入比特序列，用于后续比较
        inpData(:,bandIdx) = bitsIn; % log bits for comparison

        % Pack subband data into an OFDM symbol
        % 计算当前子带在频域中的偏移量
        offset = subbandOffset+(bandIdx-1)*subbandSize; 
        % 将子带符号填充到OFDM符号中，前后补零
        symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                         zeros(numFFT-offset-subbandSize, 1)];
        % 对填充后的OFDM符号进行IFFT变换，将其转换到时域
        ifftOut = ifft(ifftshift(symbolsInOFDM));

        % Filter for each subband is shifted in frequency
        % 为每个子带设计一个频率偏移的滤波器
        bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                     ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
        % 对IFFT输出的信号进行卷积滤波
        filterOut = conv(bandFilter,ifftOut);

        % Sum the filtered subband responses to form the aggregate transmit
        % signal
        % 将当前子带滤波后的信号累加到发射信号中
        txSig = txSig + filterOut;     
    end
    % 将当前UFMC符号的发射信号存储到复数矩阵中
    UFMCcomplex(:,symIdx) = txSig;
    % 将当前UFMC符号的所有子带输入比特序列存储到输入数据矩阵中
    inputData(:,symIdx)= inpData(:);
end

%Now we just need to convert the result to real for IM/DD
% 初始化一个实数矩阵，用于存储转换后的实数信号，列数为UFMC符号数量的两倍
UFMC_real = zeros(numFFT+filterLen-1, 2*numUFMCSymbols);
% 遍历每个UFMC符号
for symIdx = 1:numUFMCSymbols
    % 将复数信号的虚部存储到偶数列
    UFMC_real(:,(symIdx*2)) = imag(UFMCcomplex(:,symIdx));
    % 将复数信号的实部存储到奇数列
    UFMC_real(:,(symIdx*2)-1) = real(UFMCcomplex(:,symIdx));
end

%Paralel to serie and clocks
% 将实数矩阵进行并行转串行处理，得到一维的发射信号
SignalTransmitted = reshape(UFMC_real, 2*numUFMCSymbols*(numFFT+filterLen-1), 1);   


% 将输入数据保存到inputdata.mat文件中
save('inputdata.mat','inputData');
% 将发射信号保存到SignalTransmitted.mat文件中
save('SignalTransmitted.mat','SignalTransmitted');