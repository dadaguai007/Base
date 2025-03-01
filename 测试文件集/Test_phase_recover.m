%phase recover
%Transmitter parameters:
clc;clear;close all;
M   = 16 ;          % order of the modulation format
constType = 'qam';  % constellation type
Rs  = 32e9;         % symbol rate [baud]
SpS = 8 ;           % samples per symbol
Nbits = (log2(M)*1e5);   % total number of bits per polarization
pulse = 'rrc';      % pulse shaping filter
Ntaps = 1024;       % number of pulse shaping filter coefficients
alphaRRC = 0.01;    % RRC rolloff
Pch_dBm = 1;        % power per WDM channel [dBm]
Nch     = 1;        % number of WDM channels
Fc      = 193.1e12; % central optical frequency of the WDM spectrum
lw      = 100e3;    % laser linewidth
Nmodes = 2 ;        % number of signal modes [2 for polarization multiplexed signals]
freqSpac=50e9;
Tx_laser_type='ideal';
%generate WDM signal
% sigTx, symbTx_, paramTx = simpleWDMTx(paramTx);

Fs = Rs*SpS;
Ts  = 1/Rs ;
Ta  = 1/Fs;



psfRollOff=0.01;
psfLength=100;
psfShape='sqrt';


%IQM
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;
%LO
Plo_dBm  = 10;
Plo =10.^(Plo_dBm/10)*1e-3 ;
paramLO=struct();
paramLO.P = Plo;
paramLO.lw = 100e3;          % laser linewidth
paramLO.RIN_var = 0;
paramLO.Fs = Fs;
%PD
paramPD=struct();
paramPD.B =Rs;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;



% WDM each channel frequence
freqGrid = (-floor((Nch) / 2):floor((Nch-1 )/ 2)) * freqSpac;

% channel number
if mod(Nch, 2) == 0
    freqGrid = freqGrid + freqSpac / 2;
end

%power
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
        fprintf(' mode %d\t power: %.2f dBm\n', indMode, 10 * log10((Pch(indCh)/Nmodes) / 1e-3));

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
            if strcmp(Tx_laser_type,'phasenoise')
                % Generate LO field with phase noise
                phi_pn_lo = phaseNoise(lw, length(sigTx), Ta);
                sigLO = exp(1i * phi_pn_lo);
                Pin=sigLO;
            else
                Pin=10;
                Pin =10.^(Pin/10)*1e-3 ;
            end
        end
        %调制
        Amp=0.5;
        sigTxCh = iqm(Pin, Amp*sigTx,paramIQ);
        % 每个通道的发射功率是相对于所有偏振模式总功率，不是相对于单个偏振模式的功率。
        % 将发射功率除以偏振模式的数量可以得到每个通道在单个偏振模式下的功率。
        sigTxCh = sqrt(Pch(indCh) / Nmodes) * pnorm(sigTxCh);
        %mode and channel sum
        phi_sig=sigTxCh.*exp(1i*2*pi*(freqGrid(indCh)/Fs)* t);
        sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + phi_sig.';
        %total power is sum(abs(x.^2))
        Pmode = Pmode + signalpower(sigTxCh);
        %         fprintf('optical mode %d\t power: %.2f dBm\n', indMode, 10 * log10(signalpower(sigTxCh) / 1e-3));
    end

    Psig = Psig + Pmode;
    %每个通道的功率 dBm
    fprintf('channel %d\t power: %.2f dBm\n', indCh, 10 * log10(Pmode / 1e-3));
end
%整个WDM输出的功率 dBm
fprintf('total WDM signal power: %.2f dBm\n', 10 * log10(Psig / 1e-3));

% sigTxWDM
chIndex=1;

% LO
paramLO.N = length(sigTxWDM);
sigLO = basicLaserModel(paramLO);
FO      = 150e6;                 % frequency offset
f_lo   = freqGrid(chIndex)+FO;  % shift of the channel to be demodulated
sigLO = sigLO.*exp(1j*2*pi*f_lo*t); % add frequency offset

%PD
theta=0;
Sig_Rx=pdmCoherentReceiver(sigTxWDM, sigLO, theta, paramPD);
scatterplot(Sig_Rx(:,1))


%match filtering
for indMode = 1:Nmodes
Sig_Rx(:,indMode)=conv(Sig_Rx(:,indMode),hsqrt,'same');
end


% 最佳采样
paramDec = struct();
paramDec.SpS_in  = SpS;
paramDec.SpS_out = 1;
%最佳采样相位,最佳采样点为delay+1
[SigRx,delay] = decimate(Sig_Rx, paramDec);


%norm
x = pnorm(SigRx);
d = pnorm(symbTxWDM);

% 4th power frequency offset estimation/compensation
[Ei, ~] = fourthPowerFOE(SigRx, 1/Ts);
%norm
Ei = Ei./sqrt(mean(abs(Ei).^2));
for indMode = 1:Nmodes
scatterplot(Ei(:,indMode))
title('4th power frequency offset')
end

% bps

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
eyediagram(y_CPR_BPS(1:10000,2),24)
title('After BPS CPR ')


%dpll
paramCPR.alg = 'ddpll';
paramCPR.tau1 = 1/(2*pi*10e3);
paramCPR.tau2 = 1/(2*pi*10e3);
paramCPR.Kv  = 0.1;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_PLL,theta1] = cpr(x,d,paramCPR);
y_CPR_PLL = pnorm(y_CPR_PLL);
eyediagram(y_CPR_PLL(1:10000,2),24)
title('After PLL CPR ')

%viterbi
paramCPR.alg = 'viterbi';
paramCPR.N = 151;
paramCPR.Ts=Ts;
[y_CPR_viterbi,theta2] = cpr(x,d,paramCPR);
y_CPR_viterbi = pnorm(y_CPR_viterbi);
eyediagram(y_CPR_viterbi(1:10000,2),24)
title('After viterbi CPR ')
scatterplot(y_CPR_viterbi(:,indMode))
