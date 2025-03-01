% optical comb

clc;clear;close all;
SpS = 40;
Rs  = 10e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;
% 调制电压
Vpi = 2;
Vb = -Vpi/2;


Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
N=1000;
% 频率倍数因子
m=3;
f1=m*Rs;
RF1 = Creat_dither(Fs,f1,N);

%mzm
Ai= sqrt(Pi);
Amp=1;
sigTxo = mzm(Ai, Amp*RF1, Vb,Vpi);

plot_spectrum(sigTxo,Fs);

Vphi=-Vpi/2;
f2=Rs;
RF2=Creat_dither(Fs,f2,N*(f2/f1));
% outputPM = pm(sigTxo, Vphi + RF2, Vpi);
outputPM= Phase_Modulator(sigTxo, RF2, Vphi,Vpi);
plot_spectrum(outputPM,Fs);