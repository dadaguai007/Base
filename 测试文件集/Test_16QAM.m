clc;clear;close all;
addpath('D:\001-处理中\相干算法\optical_communication')
%16QAM
M=16;
SpS = 6;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W

%IQ
paramIQ=struct();
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;


%Ch
param=struct();
param.Ltotal = 40; %km
param.Lspan =20;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
param.Fs=Fs;

%linearChannel
paramCh=struct();
paramCh.Fs = param.Fs;
paramCh.L = param.Ltotal;
paramCh.D = 16;
paramCh.Fc  = 193.1e12;
paramCh.alpha = 0.2;



%PD
paramPD=struct();
paramPD.B =Rs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=Fs;


%LO  
Plo_dBm  = 10;   
f_lo = 150e6 ;  
phi_lo  = 0 ;    
lw = 500e3;
Plo =10.^(Plo_dBm/10)*1e-3 ;


% signal
bits=randi([0,1],log2(M),2^13);
% bit=bits(:);
% symbTx = modulateGray(bit, M, 'qam');
% symbTx=symbTx.';
symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso rrc típico
pulse = pulseShape('rrc', SpS, 4096, 0.01, Ts);
pulse = pulse./ max(abs(pulse));

%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(0.1,40,SpS,'sqrt');  
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');
% sigTx  = firFilter(pulse.', symbolsUp);
% save('Before modulator symbol.mat','sigTx');
% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;

figure;
plot(t(idX), real(sigTx(idX)),LineWidth=1.5)
title('TX symbol waveform')
plot_spectrum(sigTx,Fs)
title('TX symbol spectrum')

eyediagram(sigTx(1:10000),24)
title('TX symbol')
%Constellation 
scatterplot(downsample(sigTx,SpS))
title('TX symbol downsample')



%% optical modulator
Ai  = sqrt(Pi);
Amp=0.5;
sigTxo = iqm(Ai, Amp*sigTx, paramIQ);


figure;
plot(t(idX), real(sigTxo(idX)),LineWidth=1.5)
title('TXo symbol waveform')
plot_spectrum(sigTxo,Fs)
title('TXo symbol spectrum')

eyediagram(sigTxo(1:10000),24)
title('TXO symbol')
%Constellation 
scatterplot(downsample(sigTxo,SpS))
title('TXO symbol downsample')

%% oscilador local
t_lo = (0:length(sigTxo)-1) * Ta;
pn_lo  = 1*phaseNoise(lw,length(sigTxo),Ta);
sigLO = sqrt(Plo)*exp(1i*(2*pi*f_lo*t_lo + phi_lo+pn_lo));
plot_spectrum(sigLO,Fs)


%% channel 
chtype='ssfm';
if strcmp(chtype,'ssfm')
%ssfo
sigRxo=ssfm(sigTxo,param);
% sigRxo=awgn(sigRxo,20,'measured');
else
%linear
sigRxo = linearChannel(sigTxo, paramCh);
sigRxo=sigRxo.';
% sigRxo=awgn(sigRxo,20,'measured');
end

figure;hold on;
plot(t(idX), abs(sigRxo(idX)),LineWidth=1)
plot(t(idX), abs(sigTxo(idX)),LineWidth=1)
title('TXo symbol transmission')
plot_spectrum(sigTx,Fs)
title('TXo symbol spectrum transmission')
eyediagram(real(sigRxo(1:10000)),24)
title('TXo symbol transmission')
%Constellation 
scatterplot(downsample(sigRxo,SpS))
title('TXo symbol transmission downsample')


%% coherent——transmission
coherent_ipd = coherentReceiver(sigRxo,sigLO,paramPD);
coherent_ipd=pnorm(coherent_ipd);
%match filter
coherent_ipd=conv(coherent_ipd,hsqrt,'same');

%figure pd
% eyediagram(real(coherent_ipd(1:10000)),24)
% title('transmission-pd-fiber')
scatterplot(coherent_ipd)
title('transmission-pd-fiber')
scatterplot(downsample(coherent_ipd,SpS))
title('transmission-downsample-pd-fiber')

%%  coherent btb
coherent_ipd_btb = coherentReceiver(sigTxo,sigLO,paramPD);
coherent_ipd_btb=pnorm(coherent_ipd_btb);
%match filter
% coherent_ipd_btb  = firFilter(pulse.', coherent_ipd_btb);
coherent_ipd_btb=conv(coherent_ipd_btb,hsqrt,'same');

%figure pd btb
eyediagram(real(coherent_ipd_btb(1:10000)),24)
title('pd-btb')
%Constellation
scatterplot(coherent_ipd_btb)
title('btb-pd')
coherent_down_ipd_btb=downsample(coherent_ipd_btb,SpS);
scatterplot(coherent_down_ipd_btb)
title('btb-downsample-pd')



close all;

%% DSP
% CDC
paramCD = struct();
paramCD.Fs=Fs;
paramCD.L=param.Ltotal;
paramCD.D=param.D;
paramCD.Fc=param.Fc;

SigRx = cdc(coherent_ipd, paramCD);
% eyediagram(SigRx(1:10000),24)
% title('After CD compensate')
scatterplot(downsample(SigRx,SpS))
title('After CD compensate and downsample')

% downsample to two samples a symbol to equ or downsample to one symbol
paramDec = struct();
paramDec.SpS_in  = SpS;
paramDec.SpS_out = 1;
SigRx = decimate(SigRx, paramDec);

% SigRx3=downsample(SigRx,SpS);

%synchronization
[SigRx,~,ff] = sync(SigRx,symbTx.');
% scatterplot(SigRx)


% synchronization
% symbol have done the downsample from the tx  which used to DSP for reference
% symbRx = symbolSync(SigRx, symbTx.', paramDec.SpS_out, 'amp');

%norm
x = pnorm(SigRx);
d = pnorm(symbTx.');
d =d(1:length(x));

discard = 100;
% ind = 1:discard:length(SigRx);
ind = discard:length(SigRx)-discard;

% Compensate for possible phase rotation added by the channel
if 1
rot =mean(d./x);
SigRx = rot* SigRx;
scatterplot(SigRx)
%norm
x = pnorm(SigRx);
end

% % Equ
% paramEq=struct();
% paramEq.nTaps=5;
% paramEq.numIter = 1;
% % paramEq.mu =5e-3;
% paramEq.mu = [5e-3, 2e-4];
% %RLS forgetting factor
% paramEq.lambdaRLS=0.99;
% paramEq.SpS=paramDec.SpS_out;
% % coefficient taps
% paramEq.H=[];
% % Eq length
% paramEq.L = [floor(0.2*length(d)), floor(0.8*length(d))];
% % paramEq.L = [750,800];
% paramEq.Hiter=[];
% paramEq.storeCoeff='true';
% % alg is the cell
% % paramEq.alg='da-rde';
% paramEq.constSymb=pnorm(GrayMapping(M, 'qam'));
% paramEq.alg = {'da-rde','rde'};
% [yEq, H, errSq, Hiter] = mimoAdaptEqualizer(x, d, paramEq);
% scatterplot(yEq)

% 4th power frequency offset estimation/compensation
[Ei, ~] = fourthPowerFOE(SigRx, 1/Ts);
%norm
Ei = Ei./sqrt(mean(abs(Ei).^2));
scatterplot(Ei)
title('4th power frequency offset')


% Constellation
z = (0:M-1)';
y = qammod(z,M);
% CPR frequence offest
paramCPR = struct();
paramCPR.alg = 'bps';
paramCPR.N= 85;
paramCPR.B= 64;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb=y;
% paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_BPS,theta] = cpr(x,d,paramCPR);
y_CPR_BPS = pnorm(y_CPR_BPS);
eyediagram(y_CPR_BPS(1:10000),24)
title('After BPS CPR ')
if 1
scatterplot(y_CPR_BPS)
title('After BPS CPR and downsample')
else
scatterplot(downsample(y_CPR_BPS,2))
title('After BPS CPR and downsample')
end
close all;

% CPR frequence offest
paramCPR.alg = 'ddpll';
paramCPR.tau1 = 1/(2*pi*10e3);
paramCPR.tau2 = 1/(2*pi*10e3);
paramCPR.Kv  = 0.1;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_PLL,theta1] = cpr(x,d,paramCPR);
y_CPR_PLL = pnorm(y_CPR_PLL);
eyediagram(y_CPR_PLL(1:10000),24)
title('After PLL CPR ')
if 1
scatterplot(y_CPR_BPS)
title('After PLL CPR and downsample')
else
scatterplot(downsample(y_CPR_BPS,2))
title('After PLL CPR and downsample')
end

