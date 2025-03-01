clc;clear;
M=16;
SpS = 6;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;


% signal
bits=randi([0,1],log2(M),2000);
% bit=bits(:);
% symbTx = modulateGray(bit, M, 'qam');
% symbTx=symbTx.';
symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
symbTx = pnorm(symbTx);


SNR  = (-2:1:35);

% GMI  = zeros(length(SNR),1);
for i=1:length(SNR)
    snr=SNR(i);
power_x = mean(abs(symbTx).^2);
SNR_power = 10.^(snr/10) ;
sigma2=power_x/SNR_power;
noise = gaussianComplexNoise(size(symbTx), sigma2);
symbTx_awgn=symbTx+noise;
symRx=awgn(symbTx,snr,'measured');
% [GMI(i), NGMI(i)] = MonteCarlo_GMI(symbTx_awgn, symbTx, M, 'qam');
MI(i) = MonteCarlo_MI(symRx, symbTx, M, 'qam');
end

% C = log2(1 + 10.^(SNR/10));