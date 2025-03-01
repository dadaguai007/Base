s = rng(211);            % Set RNG state for repeatability
numFFT = 1024;           % Number of FFT points
numGuards = 212;         % Guard bands on both sides
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
numSymbols = 100;        % Simulation length in symbols     总符号数
bitsPerSubCarrier = 2;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 12;              % SNR in dB

% Prototype filter
%选择不同的滤波器系数HkOneSided；用于构建（FBMC）系统中的原型滤波器
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end
% Build symmetric filter
% 对称滤波器Hk。对称滤波器是FBMC系统中的关键，它有助于减少子载波之间的干扰。
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% Transmit-end processing
%   Initialize arrays

% 从FFT点数中减去两倍的保护带数量来实现。每个符号中所占用的载波数
L = numFFT-2*numGuards;  % Number of complex symbols per OFDM symbol
KF = K*numFFT;
KL = K*L;
%每个子载波的数据
dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);

sumFBMCSpec = zeros(KF*2, 1);
sumOFDMSpec = zeros(numFFT*2, 1);

% 每个符号的比特数，考虑到每个子载波的比特数和每个OFDM符号中的复数符号数量，并除以2以考虑过采样
numBits = bitsPerSubCarrier*L/2;    % account for oversampling by 2
inpData = zeros(numBits, numSymbols);
rxBits = zeros(numBits, numSymbols);
txSigAll = complex(zeros(KF, numSymbols));
symBuf = complex(zeros(2*KF, 1));
% 由于符号之间有重叠，我们需要一个缓冲区来存储当前符号及其相邻符号的重叠部分。
% 因此，缓冲区的大小至少应该是一个符号的长度加上下一个符号的长度，即 2*KF。

% Loop over symbols
for symIdx = 1:numSymbols
    
    % Generate mapped symbol data
    inpData(:, symIdx) = randi([0 1], numBits, 1);
    modData = qammod(inpData(:, symIdx), 2^bitsPerSubCarrier, ...
        'InputType', 'Bit', 'UnitAveragePower', true);
%     OQAM调制器  对于奇数符号，实部分配到奇数位置，虚部分配到偶数位置；对于偶数符号，反之。
    % OQAM Modulator: alternate real and imaginary parts
    if rem(symIdx,2)==1     % Odd symbols
        dataSubCar(1:2:L) = real(modData);
        dataSubCar(2:2:L) = 1i*imag(modData);
    else                    % Even symbols
        dataSubCar(1:2:L) = 1i*imag(modData);
        dataSubCar(2:2:L) = real(modData);
    end
    
    % Upsample by K, pad with guards, and filter with the prototype filter
    dataSubCarUp(1:K:end) = dataSubCar;
    dataBitsUpPad = [zeros(numGuards*K,1); dataSubCarUp; zeros(numGuards*K,1)]; % 填0

    X1 = filter(Hk, 1, dataBitsUpPad); %滤波
    % Remove 1/2 filter length delay
    X = [X1(K:end); zeros(K-1,1)];% 去除延时
    
    % Compute IFFT of length KF for the transmitted symbol
    txSymb = fftshift(ifft(X));
    
    % Transmitted signal is a sum of the delayed real, imag symbols
    symBuf = [symBuf(numFFT/2+1:end); complex(zeros(numFFT/2,1))]; %循环前缀
    symBuf(KF+(1:KF)) = symBuf(KF+(1:KF)) + txSymb;
    
    % Compute power spectral density (PSD)
    currSym = complex(symBuf(1:KF));
    [specFBMC, fFBMC] = periodogram(currSym, hann(KF, 'periodic'), KF*2, 1);
    sumFBMCSpec = sumFBMCSpec + specFBMC;
    
    % Store transmitted signals for all symbols
    txSigAll(:,symIdx) = currSym;
end

% Plot power spectral density
sumFBMCSpec = sumFBMCSpec/mean(sumFBMCSpec(1+K+2*numGuards*K:end-2*numGuards*K-K));
plot(fFBMC-0.5,10*log10(sumFBMCSpec));
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['FBMC, K = ' num2str(K) ' overlapped symbols'])
set(gcf, 'Position', figposition([15 50 30 30]));

BER = comm.ErrorRate;

% Process symbol-wise
for symIdx = 1:numSymbols
    rxSig = txSigAll(:, symIdx);
    
    % Add WGN
    rxNsig = awgn(rxSig, snrdB, 'measured');
    
    % Perform FFT
    rxf = fft(fftshift(rxNsig));
    
    % Matched filtering with prototype filter
    rxfmf = filter(Hk, 1, rxf);
    % Remove K-1 delay elements
    rxfmf = [rxfmf(K:end); zeros(K-1,1)];
    % Remove guards
    rxfmfg = rxfmf(numGuards*K+1:end-numGuards*K);
    
    % OQAM post-processing
    %  Downsample by 2K, extract real and imaginary parts
    if rem(symIdx, 2)
        % Imaginary part is K samples after real one
        r1 = real(rxfmfg(1:2*K:end));
        r2 = imag(rxfmfg(K+1:2*K:end));
        rcomb = complex(r1, r2);
    else
        % Real part is K samples after imaginary one
        r1 = imag(rxfmfg(1:2*K:end));
        r2 = real(rxfmfg(K+1:2*K:end));
        rcomb = complex(r2, r1);
    end
    %  Normalize by the upsampling factor
    rcomb = (1/K)*rcomb;
    
    % De-mapper: Perform hard decision
    rxBits(:, symIdx) = qamdemod(rcomb, 2^bitsPerSubCarrier, ...
        'OutputType', 'bit', 'UnitAveragePower', true);
end

% Measure BER with appropriate delay
BER.ReceiveDelay = bitsPerSubCarrier*KL;
ber = BER(inpData(:), rxBits(:));

% Display Bit error
disp(['FBMC Reception for K = ' num2str(K) ', BER = ' num2str(ber(1)) ...
    ' at SNR = ' num2str(snrdB) ' dB'])

