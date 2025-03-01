function [sigTxWDM, symbTxWDM, param] = simpleWDMTx(param)
% Implement a simple WDM transmitter.

% Check input parameters
if nargin <1
param.M = 16;
param.constType =  'qam';
param.Rs =  32e9;
param.SpS =  16;
param.Nbits = 60000;
param.pulse =  'rrc';
param.Ntaps =  4096;
param.alphaRRC =  0.01;
param.Pch_dBm =  -3;
param.Nch =5;
param.Fc = 193.1e12;
param.lw =  0;
param.freqSpac = 50e9;
param.Nmodes =  1;

end



% Transmitter parameters
Ts = 1 / param.Rs;  % Symbol period [s]
Fs = 1 / (Ts / param.SpS);  % Sampling frequency [samples/s]

% Central frequencies of the WDM channels
freqGrid = (-floor((param.Nch) / 2):floor((param.Nch-1 )/ 2)) * param.freqSpac;

if mod(param.Nch, 2) == 0
    freqGrid = freqGrid + param.freqSpac / 2;
end

if iscell(param.Pch_dBm)
    if length(param.Pch_dBm) == param.Nch
        Pch = (10.^ (cell2mat(param.Pch_dBm) / 10)) * 1e-3;
        % Optical signal power per WDM channel
    end
else
    % dBm to W
    Pch = (10 .^ (param.Pch_dBm / 10)) * 1e-3;
    Pch = Pch * ones(param.Nch,1);
end

% Time array
t = (0:((param.Nbits / log2(param.M)) * param.SpS - 1));

% Allocate arrays
sigTxWDM = complex(zeros(length(t), param.Nmodes));
symbTxWDM = complex(zeros(length(t) / param.SpS, param.Nmodes, param.Nch));
% total power  of WDM signal，初始值为0。
Psig = 0;

% Pulse shaping filter
if strcmp(param.pulse, 'nrz')
    pulse = pulseShape('nrz', param.SpS);
elseif strcmp(param.pulse, 'rrc')
    pulse = pulseShape('rrc', param.SpS, param.Ntaps, param.alphaRRC, Ts);
end
pulse = pulse / max(abs(pulse));





% 开始循环
for indCh = 1:param.Nch
    %通道数
    fprintf('channel %d\t fc : %3.4f THz\n', indCh, (param.Fc + freqGrid(indCh)) / 1e12);

    Pmode = 0;
    for indMode = 1:param.Nmodes
        %模式数
        fprintf('  mode #%d\t power: %.2f dBm\n', indMode, 10 * log10((Pch(indCh) / param.Nmodes) / 1e-3));

        % Generate random bits
        bitsTx = randi([0, 1], 1, param.Nbits);

        % Map bits to constellation symbols
        %             symbTx = modulateGray(bitsTx, param.M, param.constType);
        %You need to manually change the parameters yourself
        symbTx=qammod(bitsTx,param.M,'InputType','bit','UnitAveragePower',1);

        % Normalize symbols energy to 1
        %             symbTx = symbTx / sqrt(Es);

        symbTxWDM(:, indMode, indCh) = symbTx;

        % Upsampling
        symbolsUp = upsample(symbTx, param.SpS);

        % Pulse shaping
        sigTx = firFilter(pulse, symbolsUp);
        % Optical modulation
        if indMode == 1
            % Generate LO field with phase noise
            phi_pn_lo = phaseNoise(param.lw, length(sigTx), 1 / Fs);
            sigLO = exp(1i * phi_pn_lo);
        end
        sigTxCh = iqm(sigLO, 0.5 * sigTx);
        sigTxCh = sqrt(Pch(indCh) / param.Nmodes) * (sigTxCh/sqrt(mean(abs(sigTxCh).^2)));
        %mode and channel sum
        sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + sigTxCh.*exp(1i*2*pi*(freqGrid(indCh)/Fs)* t);

        Pmode = Pmode + signal_power(sigTxCh);
    end

    Psig = Psig + Pmode;
    fprintf('channel %d\t power: %.2f dBm\n', indCh, 10 * log10(Pmode / 1e-3));
end

fprintf('total WDM signal power: %.2f dBm\n', 10 * log10(Psig / 1e-3));
param.freqGrid = freqGrid;
end
