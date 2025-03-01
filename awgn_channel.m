function y_awgn=awgn_channel(x, snr)
% the awgn is uesd for the real signal ,which is intensity signal 
% 函数现在没有那么常用，重新写了两个实数AWGN和复数AWGN噪声添加
    %intensity signal
    power_x = mean(abs(x).^2);
    %dB-to-w
    SNR_power = 10.^(snr/10);
    noise_power = power_x / SNR_power;
    noise = sqrt(noise_power) * randn(size(x)); 
    %normal(0,1)
    y_awgn= x + noise ;

end
