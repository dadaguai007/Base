clc;clear;close all;

SpS = 6;
Rs  = 10e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;

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

% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;
figure;
plot(t(idX), sigTx(idX),LineWidth=1)
ylim([-1.5,1.5])
%psd
plot_spectrum(sigTx,Fs)




% 预补偿，预失真
self=mzm_predistortion(M, 1, 0.5,'digital');


%% Generate PAM symbols
% Outputs:
% - xd = symbols at symbol rate (1 x length(dataTX))
xd = self.a(bin2gray(symbols, 'pam', M) + 1).'; % 1 x length(dataTX)
figure;
plot(t(idX), xd(idX),LineWidth=1)
ylim([-1.5,1.5])


%预补偿效果一致，已经完成