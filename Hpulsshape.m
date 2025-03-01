function H = Hpulsshape(f,h,fs)
%计算信号脉冲形状的频率响应
%% Frequency response of PAM pulse shape
% fs = Rs*SpS;
%时域的响应h
delay = grpdelay(h, 1, 1);
%频域响应
H = freqz(h/abs(sum(h)), 1, f, fs)...
    .*exp(1j*2*pi*f/fs.*delay); % remove group delay
end