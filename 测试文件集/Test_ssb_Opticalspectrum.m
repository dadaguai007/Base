clc;clear;close all;
% 生成sin电信号
fb=40e9;
fb1=60e9;
Fs = fb*10;  % 采样率
N=10000;
t = 0:(1/Fs):N*(1/Fs)-1/Fs;  % 时间向量
f= Fs * (-0.5:1/N:0.5-1/N);
signal = cos(2*pi*fb*t)+cos(2*pi*fb1*t)+1j*sin(2*pi*fb*t)+1j*sin(2*pi*fb1*t) ;  % 10 MHz正弦信号
% signal = cos(2*pi*fb*t)+1j*sin(2*pi*fb*t);  % 10 MHz正弦信号

idx=1:512;

plot_spectrum(signal,Fs)
figure;
pwelch(signal)

Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
%IQ
paramIQ=struct();
paramIQ.Vpi=2;
paramIQ.VbI=-0.7*paramIQ.Vpi;
paramIQ.VbQ=-0.7*paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;
%PD
paramPD=struct();
paramPD.B =Fs/2;
paramPD.R =1;
paramPD.type = 'real';
paramPD.Fs=Fs;



% optical modulator
Ai  = sqrt(Pi);
Amp=0.5;
sigTxo = iqm(Ai, Amp*signal, paramIQ);
sigTxo=sigTxo.';
%
% plot_spectrum(sigTxo,Fs)
figure;
pwelch(sigTxo)
% optical spectrue
[YdBm, ly, fy]=Opticalspectrum(sigTxo, 1550e-9, f,0.0001);


%pd
ipd = pd(sigTxo, paramPD);
figure
pwelch(ipd)
plot_spectrum(ipd,Fs);
