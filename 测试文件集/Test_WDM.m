clc;clear;close all;

%WDM
Nbits = 60000;
Fc = 193.1e12;
% every channel sent power
Pch_dBm =  0;
Nch =1;
freqSpac = 50e9;
Nmodes =  1;
pulse='nrz';
psfRollOff=0.01;
psfLength=100;
psfShape='sqrt';
lw = 1e3;
% Central frequencies of the WDM channels

freqGrid = (-floor((Nch) / 2):floor((Nch-1 )/ 2)) * freqSpac;

if mod(Nch, 2) == 0
    freqGrid = freqGrid + freqSpac / 2;
end

%
if iscell(Pch_dBm)
    if length(Pch_dBm) == Nch
        Pch = (10.^ (cell2mat(Pch_dBm) / 10)) * 1e-3;
        % Optical signal power per WDM channel
    end
else
    % dBm to W
    Pch = (10 .^ (Pch_dBm / 10)) * 1e-3;
    Pch = Pch * ones(Nch,1);
end

%16QAM
M=16;
SpS = 16;
Rs  = 20e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;


%IQM
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;

Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W

%ssfm
param=struct();
param.Ltotal = 800; %km
param.Lspan =10;
param.hz= 0.1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
param.Fs=Fs;

%PD
paramPD.B =Rs;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;


%LO  
Plo_dBm  = 10;   
Plo =10.^(Plo_dBm/10)*1e-3 ;

rng(1);
% Time array
t = (0:((Nbits / log2(M))*SpS - 1));

%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(psfRollOff,psfLength,SpS,psfShape);  

% pulse shape
if strcmp(pulse, 'nrz')
    pulse = pulseShape('nrz', SpS);
elseif strcmp(pulse, 'rrc')
    pulse = pulseShape('rrc', SpS, 4096, psfRollOff, Ts);
end
pulse = pulse / max(abs(pulse));

% Allocate arrays
sigTxWDM = complex(zeros(length(t), Nmodes));
%symbol
symbTxWDM = complex(zeros(length(t) / SpS, Nmodes, Nch));
% total power of WDM signal，初始值为0。
Psig = 0;

for indCh = 1:Nch
    %通道数
    fprintf('channel %d\t fc : %3.4f THz\n', indCh, (Fc + freqGrid(indCh)) / 1e12);

    Pmode = 0;
    for indMode = 1:Nmodes
        %模式数,每个模式的输入功率
        fprintf('  mode %d\t power: %.2f dBm\n', indMode, 10 * log10((Pch(indCh)/Nmodes) / 1e-3));

        % Generate random bits
        bitsTx = randi([0, 1], log2(M), Nbits/log2(M));

        % Map bits to constellation symbols
        symbTx=qammod(bitsTx,M,'InputType','bit','UnitAveragePower',1);
        symbTx = pnorm(symbTx);
        % Normalize symbols energy to 1
        symbTxWDM(:, indMode, indCh) = symbTx;

        % Upsampling
        symbolsUp = upsample(symbTx, SpS);
        % Pulse shaping
        sigTx=conv(symbolsUp,hsqrt,'same');
%         sigTx = firFilter(pulse, symbolsUp);
        % Optical modulation
        if indMode == 1
            % Generate LO field with phase noise
            phi_pn_lo = phaseNoise(lw, length(sigTx), Ta);
            sigLO = exp(1i * phi_pn_lo);
        end
        %调制
        Amp=0.5;
        sigTxCh = iqm(sigLO, Amp*sigTx,paramIQ);
        % 每个通道的发射功率是相对于所有偏振模式总功率，不是相对于单个偏振模式的功率。
        % 将发射功率除以偏振模式的数量可以得到每个通道在单个偏振模式下的功率。
        sigTxCh = sqrt(Pch(indCh) / Nmodes) * pnorm(sigTxCh);
        %mode and channel sum
        phi_sig=sigTxCh.*exp(1i*2*pi*(freqGrid(indCh)/Fs)* t);
        sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + phi_sig.';
        %total power is sum(abs(x.^2))
        Pmode = Pmode + sum(abs(sigTxCh).^2);
    end

    Psig = Psig + Pmode;
    %每个通道的功率 dBm
    fprintf('channel %d\t power: %.2f dBm\n', indCh, 10 * log10(Pmode / 1e-3));
end
%整个WDM输出的功率 dBm
fprintf('total WDM signal power: %.2f dBm\n', 10 * log10(Psig / 1e-3));

plot_spectrum(sigTxWDM,Fs);

%ssfo
if Nmodes==1
    sigRxoWDM=ssfm(sigTxWDM,param);
    plot_spectrum(sigRxoWDM,Fs);
    fprintf('total ssfm signal power: %.2f dBm\n', 10 * log10(sum(abs(sigRxoWDM.^2)) / 1e-3));
  
else
    sigRxoWDM = manakovssfm(sigTxWDM, param);
    plot_spectrum(sigRxoWDM,Fs);
    fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(sum(abs(sigRxoWDM(:,1).^2)) / 1e-3));
    fprintf('total ssfm signal power model 2: %.2f dBm\n', 10 * log10(sum(abs(sigRxoWDM(:,2).^2)) / 1e-3));
end

%Receiver:
chIndex = 1 ;
%LO  
FO      = 150e6;
f_lo   = freqGrid(chIndex)+FO;
lw_lo      = 100e3;
phi_lo  = 0;  

pn_lo  = 1*phaseNoise(lw_lo,length(sigTx),Ta);
% oscilador local
t_lo = (0:length(sigTx)-1) * Ta;
sigLO = sqrt(Plo)*exp(1i*(2*pi*f_lo*t_lo + phi_lo+pn_lo));


%coherent receiver
if Nmodes>1
    %双偏振
    theta=pi/3;
    coherent_ipd=pdmCoherentReceiver(sigRxoWDM, sigLO, theta, paramPD);
else
    coherent_ipd = coherentReceiver(sigRxoWDM,sigLO,paramPD);
end
% 匹配滤波
for i=1:Nmodes
coherent_ipd(:,i)=conv(coherent_ipd(:,i),hsqrt,'same');
end
% coherent_ipd = firFilter(pulse, coherent_ipd);
% use rectangle filter also  can receiver the signal

plot_spectrum(coherent_ipd,Fs)
scatterplot(coherent_ipd(:,2))


%DSP
% CDC
paramCD = struct();
paramCD.Fs=Fs;
paramCD.L=param.Ltotal;
paramCD.D=param.D;
paramCD.Fc=param.Fc-f_lo;

SigRx = cdc(coherent_ipd, paramCD);
eyediagram(SigRx(1:10000),24)
title('After CD compensate')
scatterplot(downsample(SigRx,SpS))
title('After CD compensate and downsample')


% downsample to two samples a symbol to equ or downsample to one symbol
paramDec = struct();
paramDec.SpS_in  = SpS;
paramDec.SpS_out = 2;
%最佳采样相位
SigRx = decimate(SigRx, paramDec);

%synchronization
[SigRx,~,ff] = sync(SigRx,symbTx.');


%norm
x = pnorm(SigRx);
d = pnorm(symbTx.');
% d =d(1:length(x));



discard = 100;
% ind = 1:discard:length(SigRx);
ind = discard:length(SigRx)-discard;

% Compensate for possible phase rotation added by the channel
if 0
rot =mean(d./x);
SigRx = rot* SigRx;
scatterplot(SigRx)
%norm
x = pnorm(SigRx);
end

% Equ
paramEq=struct();
paramEq.nTaps=5;
paramEq.numIter = 1;
% paramEq.mu =5e-3;
paramEq.mu = [5e-3, 2e-4];
%RLS forgetting factor
paramEq.lambdaRLS=0.99;
paramEq.SpS=paramDec.SpS_out;
% coefficient taps
paramEq.H=[];
% Eq length
paramEq.L = [floor(0.2*length(d)), floor(0.8*length(d))];
% paramEq.L = [750,800];
paramEq.Hiter=[];
paramEq.storeCoeff='true';
% alg is the cell
% paramEq.alg='da-rde';
paramEq.constSymb=pnorm(GrayMapping(M, 'qam'));
paramEq.alg = {'da-rde','rde'};
[yEq, H, errSq, Hiter] = mimoAdaptEqualizer(x, d, paramEq);
scatterplot(yEq)

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