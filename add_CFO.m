function y_CFO=add_CFO(y,CFO,Nfft)
% To add an arbitrary frequency offset
% Input: y    = Time-domain received signal
%        dCFO = FFO (fractional CFO) + IFO (integral CFO)
%        Nfft = FFT size;

if isrow(y)
    % 输入向量已经是行向量，无需更改
    y = y;
else
    % 输入向量不是行向量，转换为行向量
    y = y';
end
nn=0:length(y)-1; 
y_CFO = y.*exp(1i*2*pi*CFO*nn/Nfft);