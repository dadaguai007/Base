clc;clear;close all;

SpS = 6;
Rs  = 60e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;

%MZM
Vpi = 2;
Vb = -Vpi/2;
Pi = 10^(Pi_dBm/10)*1e-3; %W

%fiber
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
% param.amp='none';
param.Fs=Fs;

%PD
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
% hsqrt=hsqrt./max(abs(hsqrt));
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



%mzm
Ai     = sqrt(Pi);
% 放大倍数
Amp=1;
sigTxo = mzm(Ai, Amp*sigTx, Vb,Vpi);


chtype='ssfm';
if strcmp(chtype,'ssfm')
%ssfo
sigRxo=ssfm(sigTxo,param);
% sigRxo=awgn(sigRxo,20,'measured');
else
%linear
sigRxo = linearChannel(sigTxo, paramCh);
sigRxo=sigRxo.';

end

%pd
ipd = pd(sigRxo, paramPD);
Ipd=conv(ipd,hsqrt,'same');
figure
plot(t(idX),Ipd(idX),LineWidth=1)
eyediagram(Ipd,2*SpS)
title('pd-after-fiber')
plot_spectrum(Ipd,Fs);

% 功率选择性衰落
[powershift,fshift,p1,fnew,p2,f] = FFT(Ipd,Fs);
figure;
plot(fnew,20*log(p1))
% 频谱表示
F=abs(fftshift(fft(Ipd)));
figure;
plot(f,20*log10(F))

%% MI
SNRs=1:30;
chi=[-3,-1,1,3];
chi=pnorm(chi);
X=symbTx;
SNR_power = 10.^(SNRs/10);
for idx=1:length(SNRs)
    Y = awgn_channel(X, SNRs(idx));
    % Y=awgn(X, SNRs(idx),'measured');
    pxy = Calc_px_given_y(X, Y, mean(abs(X).^2)/SNR_power(idx), chi);
    mi(idx) = log2(M)-(-mean(log2(pxy)));
end
awgn_Captacal=0.5 * log2(1+ SNR_power);
figure;
plot(SNRs,mi)
hold on
plot(SNRs,awgn_Captacal)