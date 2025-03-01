function [powershift,fshift,p1,fnew,p2,f] = FFT(OUT,fs)
m = length(OUT);
X = fft(OUT);
% f = (0:m-1)*(fs/m);     %frequency range
% power = abs(X).^2/m;    %power
% X = X/m;
%频谱搬移
Y = fftshift(X);
fshift = (-m/2:m/2-1)*(fs/m); % zero-centered frequency range
powershift = abs(Y).^2/m;     % zero-centered power
%更改过的频谱
p = abs(fft(OUT)/m);
p1 = p(1:m/2);
p1(2:end-1) = 2*p1(2:end-1);
fnew = (0:(m/2-1))*fs/m;
% frequency_samples(:,i)= P1; 提取频率到新的数组中的语句
%绘制更清晰的功率谱
p2=fftshift(10*log10(abs(fft(OUT))));
f= fs * linspace( -0.5, 0.5, m ).' / 1e9;
end

% FFT运算函数