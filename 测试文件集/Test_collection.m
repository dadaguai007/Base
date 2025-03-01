clc;clear;close all;
SpS=4;
Rs  = 30e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;
% PAM的整形
M=16;
symbols=(0:M-1);
% Mapeia bits para pulsos eletricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);

P=2;
% 出来的是符号分配的概率
probSymb = maxwellBoltzmann(P, symbTx);
Es = sum(abs(symbTx).^2 .*probSymb);
symbTx = symbTx ./ sqrt(Es);
figure;
stem(symbTx,probSymb)



% 假设 constSymb 和 probSymb 已经定义
numSymbols = 2000000;

SymbTx = randsample(symbTx, numSymbols, true, probSymb);
SymbTx = pnorm(SymbTx);
% Upsampling
symbolsUp = upsample(SymbTx, SpS);

% Pulso 
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
sigTx  = firFilter(pulse, symbolsUp);
if 1
    hsqrt = rcosdesign(0.01,256,SpS,'sqrt');
    % % pulse shaping
    sigTx=conv(symbolsUp,hsqrt,'same');
end
sigTx=awgn(sigTx,10,'measured');
figure;
histogram((sigTx))
if 0
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
param.D = 0;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
% param.amp='none';
param.Fs=Fs;



%mzm
Ai     = sqrt(Pi);
% 放大倍数
Amp=1;
sigTxo = mzm(Ai, Amp*sigTx, Vb,Vpi);

% sigRxo=ssfm(sigTxo,param);

figure;
histogram(real(sigTxo))
end

%%
SpS=4;
M=16;
symbols=[0:M-1];
% Mapeia bits para pulsos eletricos
symbTx = qammod(symbols,M,'gray');
symbTx = pnorm(symbTx);


P=-2.5;
% 出来的是符号分配的概率
probSymb = maxwellBoltzmann(P, symbTx);
Es = sum(abs(symbTx).^2 .*probSymb);
symbTx = symbTx ./ sqrt(Es);
figure;
stem3(real(symbTx), imag(symbTx), probSymb, 'BaseValue', 0,  'LineWidth', 1);

% 假设 constSymb 和 probSymb 已经定义
numSymbols = 200000;

SymbTx = randsample(symbTx, numSymbols, true, probSymb);
% Upsampling
symbolsUp = upsample(SymbTx, SpS);

% Pulso 
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
sigTx  = firFilter(pulse, symbolsUp);

sigTx=awgn(sigTx,20,'measured');

scatterplot(downsample(sigTx,SpS))
