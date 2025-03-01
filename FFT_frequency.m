function  [compensatedData1,estFreqOffset1,fftData]=FFT_frequency(offsetData,param)
% 使用频谱进行频偏估计

mOrder=param.M;
fs=param.Fs;

% input signal is N × 1 vector

% 频率分辨率
fr = 1;
% 根据频率分辨率计算FFT点数
Nfft = 2^(floor(log2(fs/fr)));
% 计算FFT，取幅值
% 注意要对输入数据做mOrder的乘方运算将所有星座上的频偏都对应到一个值，即mOrder倍的频偏上
fftData = abs(fft(offsetData.^mOrder,Nfft));

% 最大值索引
[~,maxIndex] = max(fftData);

% 将频率索引对应到-fs/2~fs/2上
if maxIndex>Nfft/2
    maxIndex = maxIndex - Nfft;
end

% 计算频偏。注意要除以mOrder，因为在计算FFT的时候，对输入数据进行了mOrder的乘方运算
estFreqOffset1 = fs/Nfft*(maxIndex-1)/mOrder;

% compensate the FreqOffset to the signal
compensatedData1 = offsetData .*exp(2*pi*1i*(-estFreqOffset1)/fs*(1:nSamp)');
end