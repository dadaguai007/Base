clc;clear;close all;

% OOK transimision
SpS = 6;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

Pi_dBm = 10;
% Pi_dBm=[-10:10];
Pi = 10^(Pi_dBm/10)*1e-3; %W

Vpi = 2;
Vb = -Vpi/2;


param=struct();
param.Ltotal = 40; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
param.Fs=Fs;


paramPD=struct();
paramPD.B =Rs;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;

%linearChannel
paramCh=struct();
paramCh.Fs = param.Fs;
paramCh.L = param.Ltotal;
paramCh.D = 16;
paramCh.Fc  = 193.1e12;
paramCh.alpha = 0.2;

paramEDFA = struct();
paramEDFA.G = paramCh.alpha*paramCh.L;    
paramEDFA.NF = 4.5  ;
paramEDFA.Fc = paramCh.Fc;
paramEDFA.Fs = Fs;
paramEDFA.type='noise';

bits = randi([0, 1], 1, 10000);
n = 0:length(bits)-1;

% Mapeia bits para pulsos elétricos
symbTx = 2 * bits - 1;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso NRZ típico
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));

% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);

% idX = 0:512;
idX = 1:1024;

sigTx  = firFilter(pulse, symbolsUp);


figure;
plot(t(idX), real(sigTx(idX)),LineWidth=1)
ylim([-1.5,1.5])
%psd
% plot_spectrum(sigTx,Fs);


%mzm
Ai     = sqrt(Pi);
Amp=1;
sigTxo = mzm(Ai, Amp*sigTx, Vb,Vpi);

figure
plot(t(idX), real(sigTxo(idX)),LineWidth=1)
ylim([-0.05,0.15])
%psd
% plot_spectrum(sigTxo,Fs);

%ssfm
% sigRxo=ssfm(sigTxo,param);

sigCh = edfa(sigTxo, paramEDFA);

chtype='lin';
if strcmp(chtype,'ssfm')
%ssfo
sigRxo=ssfm(sigTxo,param);
else
%linear
sigRxo = linearChannel(sigTxo, paramCh);
sigRxo=sigRxo.';

end


%pd
ipd = pd(sigRxo, paramPD);
figure
plot(ipd(idX),LineWidth=1)
eyediagram(ipd(1:10000),24)
title('pd-after-fiber')

%norm
I_Rx = ipd/std(ipd);

% capture samples in the middle of signaling intervals
symbRx = downsample(I_Rx,SpS);

%subtract DC level and normalize power
symbRx = symbRx - mean(symbRx);
symbRx = pnorm(symbRx);
M=2;
snr = signalpower(symbRx)/(2*signalpower(symbRx-symbTx));
EbN0 = 10*log10(snr/log2(M));
fprintf('the snr of signal power：%.2f dB\n',10 * log10(snr / 1e-3))

%demodulate symbols to bits with minimum Euclidean distance 

const = GrayMapping(M,'pam');% get PAM constellation
%%calculate the average energy per symbol of the PAM constellation
Es = signalpower(const);
fprintf('the power of pam constellation %.2f W\n',Es)

%
I1 = mean(I_Rx(bits==1));%  average value of I1
I0 = mean(I_Rx(bits==0));%  average value of I0

sigma1 = std(I_Rx(bits==1)); %standard deviation σ1 of I1
sigma0 = std(I_Rx(bits==0)); %standard deviation σ0 of I0

Id = (sigma1*I0 + sigma0*I1)/(sigma1 + sigma0); % optimal decision threshold
Q = (I1-I0)/(sigma1 + sigma0);

fprintf('the decision threshold Id = %.2f \n',Id)
% 还差一个解码部分


