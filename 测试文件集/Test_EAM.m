%Test EAM
clc;clear;close all;

SpS = 6;
Rs  = 60e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = -10;


Pi = 10^(Pi_dBm/10)*1e-3; %W

%fiber
param=struct();
param.Ltotal = 80; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
% param.amp='none';
param.Fs=Fs;

%PD
paramPD=struct();
paramPD.B =Rs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=Fs;

%linearChannel
paramCh=struct();
paramCh.Fs = param.Fs;
paramCh.L = param.Ltotal;
paramCh.D = 16;
paramCh.Fc  = 193.1e12;
paramCh.alpha = 0.2;

%PAM
M=4;
data_2bit=randi([0,1],log2(M),8000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;

% Mapeia bits para pulsos eletricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);
% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso 
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
% filter pulse
% sigTx  = firFilter(pulse, symbolsUp);

% %Pulso
hsqrt = rcosdesign(0.01,256,SpS,'sqrt');  
% % pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');

% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;
figure;
plot(t(idX), sigTx(idX),LineWidth=1)
ylim([-1.5,1.5])
title('Tx')
%psd
plot_spectrum(sigTx,Fs);



%DML
%besslf filter
% fcnorm=0.8*Rs/(Fs/2);
% order=3;
% filt=Besslf_filter(order,fcnorm,false);

% two_pol_filter
filt=two_pol_filter(Rs,Fs,false);
N=length(sigTx);
% f1=Fs*(-0.5:1/N:0.5-1/N);
[f, t] = freq_time_set(N, Fs);
H=filt.H;
figure;
plot(f,20*log10(abs(H(f/Fs))))


paramEAM.f=f/Fs;
%transient chirp
paramEAM.alpha=1;
paramEAM.H=H;
paramEAM.Fs=Fs;
%放大倍数
Amp=1;
sigTxo=eam(Pi,Amp*sigTx,paramEAM);