clc;clear;close all;

SpS = 4;
Rs  = 1e6;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;   
% Fs  = 2e9; 
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;

%OOK
% prbs
bits = prbs1(9,2044,1);
bits=double(bits)';

% Mapeia bits para pulsos elétricos
symbTx = 2 * bits - 1;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso NRZ típico
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
% filter pulse
sigTx  = firFilter(pulse, symbolsUp);
% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:511;
figure;
plot(t(idX), sigTx(idX),LineWidth=1)
ylim([-1.5,1.5])
%psd
plot_spectrum(sigTx,Fs)
% lowpass
% fType='gauss';
% B=20e9;
% N=8001;
% h = lowpassFIR(B, Fs, N, fType);
% ipd = firFilter(h, sigTx);
% plot_spectrum(ipd,Fs)

% fir
N=8001;
fc_hpf=5e9;
fp_bandpass=[1e9 5e9];
wn_hpf=fc_hpf*2/Fs;
wn_bandpass=fp_bandpass*2/Fs;
b_hpf=fir1(N-1,wn_hpf,'high');
b_bandpass=fir1(N-1,wn_bandpass,'bandpass');
% 频响
m_hpf=20*log(abs(fftshift(fft(b_hpf))))/log(10);
m_bandpass=20*log(abs(fftshift(fft(b_bandpass))))/log(10);
f= Fs * (-0.5:1/N:0.5-1/N);
figure;hold on;
plot(f,m_bandpass)
plot(f,m_hpf)
% filter
Y=firFilter(b_bandpass, sigTx);
% Y=filter(b_bandpass,1,sigTx);
plot_spectrum(Y,Fs)


%%
bits = prbs(7,10000);
bits=double(bits)';

% prbs time domain
idx=1:512;
figure;
plot(bits(idx))
% 相关
y=xcorr(bits-mean(bits));
figure;
plot(y)

% Mapeia bits para pulsos elétricos
symbTx = 2 * bits - 1;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso NRZ típico
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
% filter pulse
sigTx  = firFilter(pulse, symbolsUp);
% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;
figure;
plot(t(idX), sigTx(idX),LineWidth=1)
ylim([-1.5,1.5])
%psd
plot_spectrum(sigTx,Fs)

