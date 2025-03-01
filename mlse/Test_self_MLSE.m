%  Test  oneself MLSE PAM-4 or OOK
clc;clear;close all;
addpath("D:\PhD\Codebase\");
type='PAM4';
% 参数
sps = 2;
Rs = 510e9;
Fs = Rs*sps;
% 信号生成
M=2;
msgLen = 8000;
data_2bit=randi([0,1],log2(M),msgLen);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;
symbols=symbols.';
% Mapeia bits para pulsos eletricos
symbTx = pammod(symbols,M,0,'gray');

% Pulso
hsqrt = rcosdesign(0.01,256,sps,'sqrt');

% Upsampling
symbolsUp = upsample(symbTx,sps);
% % pulse shaping
sig=conv(symbolsUp,hsqrt,'same');

% 星座点
const = pammod([0:M-1],M);

% MLSE 回溯长度
tblen =  10;

% 信道
chanest = [0.986; 0.845; 0.237; 0.12345+0.31i];
msgFilt = filter(chanest,1,sig);




% 自带的mlse
eqSym = mlseeq(msgFilt,chanest,const,tblen,'rst',sps);
eqMsg = pamdemod(eqSym,M);

tran_sig=downsample(msgFilt,sps);
eqMsg_tran = pamdemod(tran_sig,M);


[nerrs_mlse ber_mlse] = biterr(data_2bit.', eqMsg)
[nerrs ber] = biterr(data_2bit.', eqMsg_tran)