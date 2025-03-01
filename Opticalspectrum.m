function [YdBm, ly, fy] = Opticalspectrum(E, lambda, f,resolution)
% 光谱函数
c = 299792458;      % speed of light
resm = resolution*1e-9;
%% Calculate optical spectrum
% Calculate spectrum
if size(E, 2) == 1 % 1 pol
    P = abs(fftshift(fft(E))).^2;
else % 2-pol
    P = abs(fftshift(fft(E(:,1)))).^2 + abs(fftshift(fft(E(:,2)))).^2;
end

% Convert OSA resolution to Hz for wavelength = lambda
resolutionHz = c/(lambda - resm/2) - c/(lambda + resm/2);

fs = 2*max(abs(f)); % sampling rate
%normalization
P = P/(fs*length(P)/resolutionHz);
%分辨率 resolutionHz 表示光谱仪的分辨率，即每单位频率间隔对应的波长范围。
% df 是实际频率分辨率。
% 通过将光谱分辨率与频率分辨率相除，可以得到需要的滤波窗口长度，以确保在计算光谱时不会失去重要的信息。
df = abs(f(2) - f(1));
windowLength = round(resolutionHz/df);
% odd number
if mod(windowLength, 2) == 0
    windowLength = windowLength + 1;
end
% Smoothing
Pfilt = filter(ones(1, windowLength)/windowLength, 1, P);
Pfilt = circshift(Pfilt, -(windowLength-1)/2); % remove delay due to filter
% center frequence
fc = c/lambda;
%fliplr反转向量，确保向量按递增的顺序排列
fy = fc + [fliplr(0:-resolutionHz:min(f)) resolutionHz:resolutionHz:max(f)];
% 对新生成的频率轴进行插值，以符合光谱的分辨率 (线性插值)
Y = interp1(f, Pfilt, fy-fc);
% wavelength
ly = c./fy;
%output
YdBm = 10*log10(Y/1e-3);

figure(132), box on
plot(ly*1e9, YdBm)
xlabel('Wavelength (nm)')
ylabel('Power (dBm)')
title(sprintf('Optical spectrum with resolution %.2f nm', resolution))