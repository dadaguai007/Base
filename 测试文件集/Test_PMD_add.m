%PMD add
clc;clear;close all;

%QPSK
M=4;
SpS = 6;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;


paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;

Pi = 10^(Pi_dBm/10)*1e-3; %W

param=struct();
param.Ltotal = 20; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
param.Fs=Fs;


paramPD.B =30e9;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;


%LO  
Plo_dBm  = 10;   
f_lo = 0 ;  
phi_lo  = 0 ;    
lw = 1e4;
Plo =10.^(Plo_dBm/10)*1e-3 ;
% typle = 'QAM-to-QPSK';
typle = 'QPSK';
if strcmp(typle,'QAM-to-QPSK')

    bits=randi([0,1],log2(M),8000);
    symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
    symbTx = pnorm(symbTx);
else
    bits=randi([0,M-1],1,8000);
    symbTx=pskmod(bits,M,pi/4,'gray') ;
    symbTx = pnorm(symbTx);
end

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso rrc típico
% pulse = pulseShape('rrc', SpS, 4096, 0.01, Ts);
% pulse = pulse./ max(abs(pulse));
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));

hsqrt = rcosdesign(0.01,20,SpS,'sqrt');  %设计根升余弦脉冲成型滤波器
% sigTx=conv(symbolsUp,hsqrt,'same');

sigTx  = firFilter(pulse.', symbolsUp);

% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;

figure;
plot(t(idX), real(sigTx(idX)))
%psd
plot_spectrum(sigTx,Fs);

Ai  = sqrt(Pi);
Amp=0.5;
sigTxo = iqm(Ai, Amp*sigTx, paramIQ);

% figure;
% plot(t(idX), real(sigTxo(idX)))
plot_spectrum(sigTxo,Fs);


pn_lo  = phaseNoise(lw,length(sigTxo),Ta);
% oscilador local
t_lo = (0:length(sigTxo)-1) * Ta;
sigLO = sqrt(Plo)*exp(1i*(2*pi*f_lo*t_lo + phi_lo+pn_lo));



% coherent——Receiver
coherent_ipd = coherentReceiver(sigTxo,sigLO,paramPD);
coherent_ipd=pnorm(coherent_ipd);
% coherent_ipd  = firFilter(pulse.', coherent_ipd);
% coherent_ipd=conv(coherent_ipd,hsqrt,'same');

plot_spectrum(coherent_ipd,Fs);
scatterplot(coherent_ipd)

coherent_down_ipd=downsample(coherent_ipd,SpS);
scatterplot(coherent_down_ipd)



% Einf = fftshift(fft(sigTxo));
Nfft = length(sigTxo);
paramPMD.L=10e3;
paramPMD.meanDGDps=0.1;
paramPMD.PMD_span_length = 1e3;
f = Fs * (-0.5:1/Nfft:0.5-1/Nfft);
[M,N] = PMDJonesMatrix(2*pi*f,paramPMD);
% sigTxo=sigTxo*M(:,:,1);
Einf = fftshift(fft(sigTxo));
Einf = [Einf; zeros(size(Einf))];
% Ein 应为 2*N
%对列进行fft
% Einf = fftshift(fft(sigTxo));
% Nfft = length(sigTxo);
% f = Fs * (-0.5:1/Nfft:0.5-1/Nfft);
%PMD
for k = 2:length(f)
    Einf(:, k) = M(:,:,k)*Einf(:, k);
end
%对行求ifft
Eout = ifft(ifftshift(Einf,2),2);
close all;
scatterplot(Eout(1,:))
scatterplot(sigTx)
%PDL
% if size(Ein,1)==2 && PDL ~= 0
%     a = 10^(-PDL/10);
%     Eout(2, :) = a*Eout(2, :);
% end
