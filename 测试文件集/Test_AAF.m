% Test AAF
% 抗混叠滤波器
clc;clear;close all;
% 适用于 Rs（波特率）>> Fs 的情形，减小信号之间的混叠，提高信号质量
% 参数
sps = 10;
Rs = 510e9;
Fs = 256e9;
% 信号生成
M=4;
data_2bit=randi([0,1],log2(M),8000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;
symbols=symbols.';
% Mapeia bits para pulsos eletricos
symbTx = pammod(symbols,M,0,'gray');
% symbTx = pnorm(symbTx);
% Upsampling
symbolsUp = upsample(symbTx,sps);

% Pulso 
hsqrt = rcosdesign(0.01,256,sps,'sqrt');  
% % pulse shaping
sig=conv(symbolsUp,hsqrt,'same');

% AAF
aaf_sig = AAF(sig, Fs, Rs*sps);

% Resampling_Out
rs_sig = resample(aaf_sig, Fs, Rs*sps);


figure;hold on;
plot(real(aaf_sig))
plot(real(rs_sig))