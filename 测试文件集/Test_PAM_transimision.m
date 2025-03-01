clc;clear;close all;
% PAM transimision

SpS = 6;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

Vpi = 2;
Vb = -Vpi/2;


Pi_dBm = 10;
% Pi_dBm=[-10:10];
Pi = 10^(Pi_dBm/10)*1e-3; %W


param=struct();
% param.Ltotal=[0:20];
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

%PAM
M=4;
data_2bit=randi([0,1],log2(M),80000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;

% Mapeia bits para pulsos elétricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);
% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso 
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
% filter pulse
sigTx  = firFilter(pulse, symbolsUp);

% %Pulso
% hsqrt = rcosdesign(0.01,256,SpS,'sqrt');  
% hsqrt=hsqrt./max(abs(hsqrt));
% % pulse shaping
% sigTx=conv(symbolsUp,hsqrt,'same');


%mzm
Ai= sqrt(Pi);
Amp=0.25;
sigTxo = mzm(Ai, Amp*sigTx, Vb,Vpi);
power=signalpower(sigTxo);
% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;
figure;
plot(t(idX), sigTx(idX),LineWidth=1)
ylim([-1.5,1.5])
%psd
plot_spectrum(sigTx,Fs);

% channal
sigCh = edfa(sigTxo, paramEDFA);
sigRxo = linearChannel(sigCh, paramCh);
sigRxo=sigRxo.';

%pd
ipd = pd(sigRxo, paramPD);
figure;
plot(t(idX), ipd(idX),LineWidth=1)


%norm
I_Rx = ipd/std(ipd);

% capture samples in the middle of signaling intervals
symbRx = downsample(I_Rx,SpS);

%subtract DC level and normalize power
symbRx = symbRx - mean(symbRx);
symbRx = pnorm(symbRx);

snr = signalpower(symbRx)/(2*signalpower(symbRx-symbTx));
EbN0 = 10*log10(snr/log2(M));
fprintf('the snr of signal power：%.2f dB\n',10 * log10(snr / 1e-3))

%demodulate symbols to bits with minimum Euclidean distance 
const = GrayMapping(M,'pam');% get PAM constellation
%%calculate the average energy per symbol of the PAM constellation
Es = signalpower(const);
fprintf('the power of pam constellation %.2f mW\n',Es)

% 还差一个解码部分