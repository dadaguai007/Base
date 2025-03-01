% Test DAC and filter_ele
clc;clear;close all;

% 参数
sps = 10;
Rs = 510e9;
Fs = Rs*sps;
bwl = Rs*0.44;
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

% DAC
param.snr_dB = inf;
param.quantizer = struct('bits', 8, 'ENoB', 16, 'type', 'tread');
param.dac = struct('targetFs', Fs, 'Fs', Fs);
param.efilter = struct('type', 'gauss', 'bwl', bwl, 'order', 5,'sps',sps,'outputVoltage',0);
drive = DAC(sig, param);


figure;hold on;
plot(sig)
plot(drive)