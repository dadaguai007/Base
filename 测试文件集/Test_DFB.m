clc;clear ;
Len = 2^13;
Baudrate = 0.622e9;
Up = 10;
fs = Baudrate * Up;
m=0.8;
rng(2)
TxSym = randi([0 1],1,Len);
TxWfm = reshape(repmat(TxSym,Up,1),1,[]);
symLen=2^14;
TxWfm=TxWfm(1:symLen);
 
%% DML parameter
alpha = 2;
k = 1e12;
Pmax = 1e-3;  % Unit: W
 
 
%% DML modulation 
d = 1-m+m*TxWfm;  % modified according to the paper
Phase = alpha/2*log(d) + alpha/2*k*Pmax*cumsum(d)*1/fs;
Eout = sqrt(Pmax.*d).*exp(1i*Phase);
 
 
plot(linspace(-fs/2/1e9, fs/2/1e9, symLen), 10*log10(abs(fftshift(fft(Eout))).^2 /symLen));
xlabel('Frequency/GHz');ylabel('Power/dB');
title('matlab')
grid on;