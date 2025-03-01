addpath('Functions')  
clc;clear;
msg = randi([0 1],1000,1);
symb = nrModuMapper(msg,'16QAM');
% M=16;
% bits=randi([0,1],log2(M),8000);
% symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
% msg=bits(:);
N0 = 0.1;
% rxsymb = symb + sqrt(N0/2)*randn(size(symb)) ;
rxsymb = awgn(symb,1/N0,1,'linear');
msg_hat = nrSoftModuDemapper(rxsymb,'16QAM',N0,'max-log-map');
numErr = sum(msg ~=(msg_hat < 0))
plotNrConstellation('16QAM')
plotNrLLR('16QAM')