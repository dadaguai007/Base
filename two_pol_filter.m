function   filt=two_pol_filter(bw,fs,verbose)
%归一化的截止频率
% verbose=1;
type='two-pole';
% Variables used when calculating impulse response
maxMemoryLength = 2^10; % maximum memory length
threshold = 1-1e-6;      % energy threshold above which impulse respone is neglected

% 2nd-order filter with unit damping
% The continuous-time transfer function is converted to
% discrete-time using the bilinear transformation with frequency
% prewarping

% bw = 10e9; % the argument order is used as bandwidth
% fs = 50e9; % the argument fcnorm is used as sampling frequency


fcnorm = bw/(fs/2);
order = 2;


wc = 2*pi*bw/sqrt(sqrt(2)-1); % converts BW into fc
bs = 1;
as = [(1/wc).^2, 2/wc, 1];
[num, den] = bilinear(bs, as, fs, bw);

if verbose
    %     verbose = false;
    f = linspace(0, fs/2);
    % 模拟滤波器的响应
    Hs = freqs(bs, as, 2*pi*f);
    % 数字滤波器的响应
    Hz = freqz(num, den, f, fs);
    plot_transform_num_analog(Hs, Hz, f/fs, bw/fs);
end
filt.type = type;
filt.order = order;
filt.num = num;
filt.den = den;
filt.grpdelay = grpdelay(num, den, 1);
filt.fcnorm = fcnorm;
filt.H = @(f) freqz(num, den, 2*pi*f).*exp(1j*2*pi*f*filt.grpdelay);

if den == 1 % FIR
    filt.h = num;
else
    x = zeros(1, maxMemoryLength+1);
    x(1) = 1;
    y = filter(filt.num, filt.den, x);
    E = cumsum(abs(y).^2)/sum(abs(y).^2);
    y(E > threshold) = [];
    y = y/abs(sum(y)); % normalize to have unit gain at DC
    filt.h = y;
end
end