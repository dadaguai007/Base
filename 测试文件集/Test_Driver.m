% Test Amp Driver
clc;clear;
addpath('D:\PhD\Project\单边带光发射机自适应偏压控制\Simulation\传输_影响')
OFDM_TX;

[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
xt=signal;
Vgain=4;
Vbiasadj=1;
out=Driver(xt,Vgain,Vbiasadj);

figure;
subplot(2,1,1)
plot(real(signal))
subplot(2,1,2)
plot(real(out))