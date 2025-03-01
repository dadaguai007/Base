% 频域均衡器等操作（待完善）。。。
function Nops = freq_domain_equalizer(Nfft, M)
fft_type = 'split radix';
switch lower(fft_type)
    case 'split radix'
        Offt = @(Nfft) 4*Nfft.*log2(Nfft); % number of real operations assuming split-radix FFT
    otherwise
        error('Nops_freq_domain_equalizer: invalid fft_type')
end
% Number of operations divided by number of useful samples
Nops = (2*Offt(Nfft) + Nfft)./(Nfft - M);

%计算每个有用样本的总操作数。
% 2*Offt(Nfft) 表示分裂基数 FFT 中两个操作的总操作数（正向和反向）。
% Nfft 表示乘法操作的数量（在频域中的复共轭乘法）。
% (Nfft - M) 是有用样本的数量（FIR 滤波器的输出）。
% 结果：总操作数与有用样本数的比率。

end