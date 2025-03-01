% Boss 张系统测试
clc;clear;close all;

%QPSK
M=4;
SpS = 6;
Rs  = 10e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;

% 偏振数
Nmodes = 2 ;

% pulse
psfRollOff=0.01;
psfLength=100;
psfShape='sqrt';

%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(psfRollOff,psfLength,SpS,psfShape);

% 调制器设定
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;
Pi_dBm = 10;
Pin = 10^(Pi_dBm/10)*1e-3; %W

% PD
paramPD.B =30e9;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;


% 偏振调制

for indMode = 1:Nmodes
    % 调制
    rng(10);
    % typle = 'QAM-to-QPSK';
    typle = 'QPSK';
    if strcmp(typle,'QAM-to-QPSK')

        bits=randi([0,1],log2(M),8000);
        symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
        symbTx = pnorm(symbTx);
    else
        bits=randi([0,M-1],1,8000);
        symbTx=pskmod(bits,M,pi/4,'gray') ;
        symbTx = pnorm(symbTx);
    end
    % Normalize symbols energy to 1
    symbTx_Dual(:, indMode) = symbTx;
    % Upsampling
    symbolsUp = upsample(symbTx, SpS);
    % Pulse shaping
    sigTx=conv(symbolsUp,hsqrt,'same');


    % 双偏调制
    Amp=0.5;
    if indMode==1
        sigTxCh = iqm(Pin, Amp*sigTx,paramIQ);
    else
        sigTxCh = DP_iqm(Pin, Amp*sigTx,paramIQ);
    end
    sigTx_Dual(:, indMode) =sigTxCh;
end


% 接收

% 时间戳
[~,t_up]=freq_time_set(8000*SpS,Fs);

% LO
Plo_dBm  = 10;
Plo =10.^(Plo_dBm/10)*1e-3 ;
paramLO=struct();
paramLO.P = Plo;
paramLO.lw = 100e3;          % laser linewidth(相噪的产生)
paramLO.RIN_var = 0;
paramLO.Fs = Fs;
paramLO.N = length(sigTxCh);

sigLO = basicLaserModel(paramLO);
% frequency offset
f_lo   = 0;  
sigLO = sigLO.*exp(1j*2*pi*f_lo*t_up); % add frequency offset

% 双偏PD接收
theta=0;
Sig_Rx=pdmCoherentReceiver(sigTx_Dual, sigLO, theta, paramPD);

Sig_Rx_Dual=ones(size(Sig_Rx));
%match filtering
for indMode = 1:Nmodes
Sig_Rx_Dual(:,indMode)=conv(Sig_Rx(:,indMode),hsqrt,'same');
Sig_Rx_Dual(:,indMode)=pnorm(Sig_Rx_Dual(:,indMode));
end
