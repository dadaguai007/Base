function freq_offset = fft_freq_offset(sig, Fs, average_over_modes, fft_size)
    % Find the frequency offset by searching in the spectrum of the signal
    % raised to 4. Doing so eliminates the modulation for QPSK but the method also
    % works for higher order M-QAM.

    % Adjust fft_size to be a power of 2

    if (~mod(log2(fft_size),2) == 0) || (mod(log2(fft_size),2) == 1 )
        fft_size = 2^(ceil(log2(fft_size)));
    end

    % Fix number of stuff

    [npols, L] = size(sig);

    % Find offset for all modes
    freq_sig = zeros(npols, fft_size);
    for l = 1:npols
        freq_sig(l, :) = abs(fft(sig(l, :).^4, fft_size)).^2;
    end

    % Extract corresponding FO
    freq_offset = zeros(npols, 1);
    freq_vector=Fs * (-0.5:1/fft_size:0.5-1/fft_size)/4;
%     freq_vector = fftfreq(fft_size, 1/os)/4;
    for k = 1:npols
        [~, max_freq_bin] = max(abs(freq_sig(k, :)));
        freq_offset(k, 1) = freq_vector(max_freq_bin);
    end

    if strump(average_over_modes,'true')
        freq_offset = mean(freq_offset) * ones(size(freq_offset));
    end
end

%
% 这段代码是一个MATLAB函数，名为`fft_freq_offset`，用于通过分析信号的频谱来估计信号的频率偏移。这种频率偏移的估计方法适用于QPSK和更高阶的M-QAM调制信号。
% 函数`fft_freq_offset`接受四个输入参数：
% - `sig`：输入信号，是一个二维数组，其中每一行代表一个不同的模式或极化状态。
% - `Fs`：信号的采样率。
% - `average_over_modes`：一个逻辑变量，如果为`true`，则函数会计算所有模式的平均频率偏移；如果为`false`，则分别计算每个模式的频率偏移。
% - `fft_size`：进行快速傅里叶变换（FFT）的大小。
% 函数的返回值`freq_offset`是一个数组，包含每个模式的频率偏移估计。
% 函数的工作流程如下：
% 1. 调整`fft_size`，确保其为2的幂次方，这样可以提高FFT的效率。如果`fft_size`不是2的幂次方，则将其调整为大于等于`fft_size`的最小2的幂次方。
% 2. 确定信号的大小`[npols, L]`，其中`npols`是模式或极化状态的数量，`L`是每个模式的样本数。
% 3. 初始化一个用于存储FFT结果的数组`freq_sig`。
% 4. 对每个模式执行以下操作：
%    - 计算`sig`的第四次幂，然后对其进行FFT，并取其绝对值的平方，结果存储在`freq_sig`中。
% 5. 初始化一个用于存储频率偏移估计的数组`freq_offset`。
% 6. 创建一个频率向量`freq_vector`，它包含了FFT结果中每个 bin 对应的频率。
% 7. 对每个模式执行以下操作：
%    - 找到`freq_sig`中能量最高的频谱 bin，即最大值的位置。
%    - 将对应的频率值存储在`freq_offset`中。
% 8. 如果`average_over_modes`为`true`，则计算所有模式的平均频率偏移，并将这个值赋给所有模式。
% 9. 返回频率偏移估计`freq_offset`。
% 总的来说，这个函数通过分析信号频谱中的最大值来估计信号的频率偏移，这种方法对于QPSK和M-QAM调制信号都有效。通过计算信号的四次幂，可以消除QPSK调制的影响，从而更容易地识别频率偏移。
