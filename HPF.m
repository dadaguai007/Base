function dataout=HPF(datain, samplerate, HighF)
    N = numel(datain);
    data = reshape(datain, 1, numel(datain));
    t = 1:1:N;
    freq = samplerate * t / N;
    filterResponse = zeros(size(freq));
    
    % 确定截止频率位置
    cutoffPosition = find(freq >= HighF);
    
    % 构建高通滤波器响应，高于截止频率的分量通过，低于截止频率的分量阻止
    filterResponse(min(cutoffPosition):N-min(cutoffPosition)-1) = 1;
    
    % 对输入信号进行快速傅里叶变换（FFT）
    dataFFT = fft(data);
    
    % 应用滤波器响应
    dataoutFFT = dataFFT .* filterResponse;
    
    % 进行逆快速傅里叶变换（IFFT）以恢复时域信号
    dataout = ifft(dataoutFFT);
    
    % 重塑输出信号为原始输入信号的形状
    dataout = reshape(dataout, size(datain));
end