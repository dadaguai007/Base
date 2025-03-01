% 矩形抗混叠滤波（Anti-Aliasing Filter）
% 以减少重采样时的混叠效应，用于对输入信号进行滤波操作，以减少在信号重采样过程中可能产生的混叠效应
%  out = antiAliasingFilter(obj, in, Nss, type)
%方法首先根据输入参数 type 来确定是上采样还是下采样操作。
% 根据采样率 Nss 和一些常数，方法计算出滤波器的频域表示 AAF。AAF 是一个时域为正弦函数的窗口函数，通常称为“Sinc 函数”。
%
% 计算变量 shiftLength，它表示将滤波器进行循环移位（circular shift）的样本数，以避免滤波延迟。
%
% 将 AAF 的长度扩展为与输入信号 in 相同长度，并在末尾填充零，以确保足够的滤波器长度。
%
% 对滤波器 AAF 进行循环移位，以准备进行滤波操作。
%
% 使用频域卷积（FFT 和 IFFT）的方式，将滤波器 AAF 应用于输入信号 in，得到滤波后的输出信号。
% AAF 是一个时域为Sinc函数的窗口函数，通常称为"Sinc函数"。
% 长度为 (Nss/obj.resamplingRate)*202+1 或 (Nss*obj.resamplingRate)*202+1，具体取决于是上采样还是下采样操作
% 从-101到101的长度分布的一个sinc时域信号,频域上就是矩形滤波器
%resamplingRate 用于缩放Sinc函数的幅度
if strcmp(type, 'upsampling')
    % AAF means 'A'nti 'A'liasing 'F'ilter
    AAF = (resamplingRate)*sinc(Nss*linspace(-101, 101, (Nss/resamplingRate)*202+1)).';
elseif strcmp(type, 'downsampling')
    AAF = (resamplingRate)*sinc(Nss*linspace(-101, 101, (Nss*resamplingRate)*202+1)).';
end
%它表示将滤波器进行循环移位（circular shift）的样本数，以避免滤波器引入的信号延迟。
shiftLength = floor(length(AAF)/2);  % Computing number of samples for circshift()
%扩展抗混叠滤波器 AAF 的长度，使其与输入信号 in 具有相同的长度，并在末尾填充零，以确保有足够的滤波器长度。
AAF(length(in)) = 0; % Filling up the filter with zeros
%循环移位
%使用 circshift 函数将滤波器 AAF 进行循环移位，以准备进行滤波操作。这是为了消除滤波器引入的信号延迟。
AAF = circshift(AAF, [-shiftLength-1, 0]); % Shifting circularly the time-domaing filter to avoid lag
out = ifft(fft(AAF).*fft(in)); % Applying the filter in frequency domain


%note:
% 滤波器进行循环移位（circular shift）的样本数为信号的AAF长度的一半是因为这是为了避免引入信号的滤波延迟。
% 
% 抗混叠滤波器（Anti-Aliasing Filter，AAF）通常是一种具有对称性的滤波器，它的时域表示是Sinc函数。Sinc函数在时域上有一个中心峰值，因此滤波器的频谱中心也位于频谱的中心。当应用这样的滤波器时，它会引入信号的一定延迟，导致输出信号在时间上相对于输入信号发生移位。
% 
% 为了抵消这个引入的延迟，通常会在滤波器的时域表示中进行循环移位，将滤波器的中心与输入信号的中心对齐。由于Sinc函数的对称性，循环移位滤波器的样本数通常是滤波器长度的一半，以确保中心对齐。
% 
% 这个循环移位的操作有助于保持信号的时间同步性，以便在重采样或滤波过程中不引入不必要的时间延迟，从而保持信号的相对位置不变。这对于信号处理和通信系统中的精确同步非常重要。

%我们已经生成了一个时域表示为Sinc函数的 AAF。然而，由于Sinc函数通常是一个无限延伸的函数，它在时域上具有无限长度。但在实际应用中，我们需要将滤波器的长度限制为与输入信号相匹配的有限长度。

% 因此，通过 AAF(length(in)) = 0，我们将 AAF 的最后一个样本（即超出输入信号长度的部分）设置为零。这样，AAF 的长度将与输入信号的长度相匹配，而不会引入额外的延迟或混叠效应。

% 这个操作确保了滤波器在滤波操作中的正确行为，使其能够处理与输入信号相同长度的信号，并减少了对信号质量的不良影响。