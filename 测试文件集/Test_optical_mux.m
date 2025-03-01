% Test_optical_mux
clc;clear;close all;

SpS = 6;
Rs  = 60e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

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
hsqrt = rcosdesign(0.01,256,SpS,'sqrt');  
%  pulse shaping
sigTx=conv(symbolsUp.',hsqrt,'same');

xt=awgn(sigTx,20,'measured');
BW_MUX=65e9;
N=2;
[yt, Ndelay] = optical_mux(2*BW_MUX, N, Ta, xt);
plot_spectrum(xt,Fs);
plot_spectrum(yt,Fs);
plot_spectrum(sigTx,Fs);