clc;clear;close all;
addpath('D:\001-处理中\相干算法\optical_communication')
%Test QPSK matlab with VPI
M=64;
SpS = 16;
Rs  = 10e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;
Nbits = 4096;   % total number of bits per polarization
% the power is the input or output
Pch_dBm = 0;        % power per WDM channel [dBm]
% mutil channel power dBm
% Pch_dBm={-50,0,-50};
Channel_power_type='input'; % input or output
Nch     = 5;        % number of WDM channels
Fc      = 193.1e12; % central optical frequency of the WDM spectrum
lw      = 100e3;    % laser linewidth
Nmodes = 1 ;        % number of signal modes [2 for polarization multiplexed signals]
freqSpac=20e9;
Tx_laser_type='ideal';
%generate WDM signal
% sigTx, symbTx_, paramTx = simpleWDMTx(paramTx);

paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;

%Amp
paramAmp = struct();
paramAmp.G = 16.82;
paramAmp.NF = 4.5;
paramAmp.Fc = 193.1e12;
paramAmp.Fs = Fs;
paramAmp.type='none';

param=struct();
param.Ltotal = 200; %km
param.Lspan =20;
param.hz= 0.1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=2*Fs;


% bits=randi([0,1],log2(M),8000);
% symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
% symbTx = pnorm(symbTx);
% symbolsUp = upsample(symbTx, SpS);


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


% Time array
t = (0:(Nbits*SpS - 1));

hsqrt = rcosdesign(0.5,8,SpS,'normal');  %设计根升余弦脉冲成型滤波器

% 频域上采样
m=2;
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
       rng(1);
%         % Generate random  bits
        bitsTx = randi([0, 1], log2(M), Nbits);
% 
%         % Map bits to constellation symbols
        symbTx=qammod(bitsTx,M,'InputType','bit','UnitAveragePower',1);
        symbTx = pnorm(symbTx);
        % Normalize symbols energy to 1
%         symbTx=symbTx(1:4096);
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
                % this is important，理想的输入模式
                Pin=1;
%                 Pch = (10 .^ (Pch_dBm / 10)) * 1e-3;
%                 Pin=sqrt(Pch);
            end
        end
        %调制
        Amp=1.35;
        sigTxCh = iqm(Pin, Amp*sigTx,paramIQ);
        % 每个通道的发射功率是相对于所有偏振模式总功率，不是相对于单个偏振模式的功率。
        % 将发射功率除以偏振模式的数量可以得到每个通道在单个偏振模式下的功率。
        if strcmp(Channel_power_type,'output')
            % 这句的作用是，输出的光功率 就是 设置的功率，不用进行放大了，也即前面设置的功率为输出的光功率
            sigTxCh = sqrt(Pch(indCh) / Nmodes) * pnorm(sigTxCh);
        elseif strcmp(Channel_power_type,'input')
            % 这里的作用是， 前面设置的即为输入光功率为多少，
            sigTxCh = sqrt(Pch(indCh) / Nmodes)* (sigTxCh) ;
        end
        power=signalpower(sigTxCh);
        fprintf(' mode %d optical signal power: %.2f dBm\n',indMode, 10 * log10(power / 1e-3));
        %mode and channel sum
        phi_sig=sigTxCh.*exp(1i*2*pi*(freqGrid(indCh)/Fs)* t);
        %时域相加
        sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + phi_sig.';
        %频域上采样，在频域上相加
%         N=length(phi_sig);
%         X=fft(phi_sig);
%         X_up=[X(1:N/2)  zeros(1,N*m-N)  X(N/2+1:end)];
%         x_up=m*ifft((X_up),m*N);
%         sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + X_up.';
%         sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + x_up.';
        %total power is mean(abs(x.^2))
        Pmode = Pmode + signalpower(sigTxCh);
        %         fprintf('optical mode %d\t power: %.2f dBm\n', indMode, 10 * log10(signalpower(sigTxCh) / 1e-3));
    end

    Psig = Psig + Pmode;
    %每个通道的功率 dBm
    fprintf('channel %d\t power: %.2f dBm\n', indCh, 10 * log10(Pmode / 1e-3));
end
%整个WDM输出的功率 dBm

% sigTxWDM=m*ifft(sigTxWDM,m*N);
%整体上采样
N=length(sigTxWDM);
X=fft(sigTxWDM);
X_up=[X(1:N/2).'  zeros(1,N*m-N)  X(N/2+1:end).'];
x_up=m*ifft((X_up),m*N);
sigTxWDM=x_up;
fprintf('total WDM signal power: %.2f dBm\n', 10 * log10(Psig / 1e-3));

% chIndex=2;
% % f_lo should be negetival with the center of signal frequency
% f_lo=-(freqGrid(chIndex)/(2*Fs));
% t= (0:length(sigTxWDM)-1);
% Elo = exp(1i*2*pi*f_lo*t);
% sigTxWDM=sigTxWDM.*Elo;

idX=1:1024;

power=signalpower(sigTxWDM);
fprintf('WDM signal power: %.2f dBm\n', 10 * log10(power / 1e-3));
%edfa
if 1
    %10 11 12(2000) 13
        paramAmp.G = 8-10 * log10(power / 1e-3);
%         paramAmp.G = 0-10 * log10(power / 1e-3);
%     paramAmp.G = 15-10 * log10(power / 1e-3);
    E = edfa(sigTxWDM, paramAmp);
    power1=signalpower(E);
    fprintf('signal power: %.2f dBm\n', 10 * log10(power1 / 1e-3));
    sigTxo=E;
end



sigRxo=ssfm(sigTxo,param);
power2=signalpower(sigRxo);

fprintf(' after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));
% load('1200km_8dBm_WDM.mat')

% load('WDM_test.mat')
load('WDM_amp10_200.mat')
I=inputOfiber.bands;
R=I{1,1};
% plot_spectrum(R.E,2*Fs)

[~,T]=freq_time_set(Nbits*2*SpS,Fs*2);


figure;hold on;
plot(T,real(sigRxo),'LineWidth', 1)
plot(T,real(R.E),'LineWidth', 1)
xlabel('Time/s')
ylabel('Amplitude / V')
legend('matlab','vpi','FontName', 'Arial','FontSize',14,'Location', 'best')
set(gca, 'FontName', 'Arial', 'FontSize', 12);
% 设置X轴标签字体属性
set(get(gca, 'XLabel'), 'FontName', 'Arial', 'FontSize', 14);
% 设置Y轴标签字体属性
set(get(gca, 'YLabel'), 'FontName', 'Arial', 'FontSize', 14);
title('Comparison of real part of the signal')
box on;

figure;hold on;
plot(T,imag(sigRxo),'LineWidth', 1)
plot(T,imag(R.E),'LineWidth', 1)
xlabel('Time/s')
ylabel('Amplitude / V')
legend('matlab','vpi','FontName', 'Arial','FontSize',14,'Location', 'best')
set(gca, 'FontName', 'Arial', 'FontSize', 12);
% 设置X轴标签字体属性
set(get(gca, 'XLabel'), 'FontName', 'Arial', 'FontSize', 14);
% 设置Y轴标签字体属性
set(get(gca, 'YLabel'), 'FontName', 'Arial', 'FontSize', 14);
title('Comparison of imaginary part of the signal')
box on;

power3=signalpower(R.E);
fprintf('signal power: %.2f dBm\n', 10 * log10(power3 / 1e-3));

mse_real=real(sigRxo)-real(R.E);
mse_imag=imag(sigRxo)-imag(R.E);
REAL=mean(mse_real).^2;
IMAG=mean(mse_imag).^2;
fprintf('MSE: %.2f \n', REAL);
%%
plot_spectrum(sigTxWDM,2*Fs);
% dataout=conv(sigTxWDM,hsqrt,'same');
dataout=LPF(sigTxWDM,2*Fs,30e9);
% plot_spectrum(dataout,2*Fs);
load('WDM_60Ghz.mat')
I=inputOfiber.bands;
R=I{1,1};
% dataout_vpi=conv(R.E,hsqrt,'same');
dataout_vpi=LPF(R.E,2*Fs,30e9);
% plot_spectrum(dataout_vpi,2*Fs);
% dataout=sigTxWDM;
% dataout_vpi=R.E;
plot_spectrum(dataout_vpi,2*Fs);
figure;hold on;
plot(t(1:10000),real(dataout(1:10000)),'LineWidth', 1)
plot(t(1:10000),real(dataout_vpi(1:10000)),'LineWidth', 1)
xlabel('Time/s')
ylabel('Amplitude / V')
legend('matlab','vpi','FontName', 'Arial','FontSize',14,'Location', 'best')
set(gca, 'FontName', 'Arial', 'FontSize', 12);
% 设置X轴标签字体属性
set(get(gca, 'XLabel'), 'FontName', 'Arial', 'FontSize', 14);
% 设置Y轴标签字体属性
set(get(gca, 'YLabel'), 'FontName', 'Arial', 'FontSize', 14);
title('Comparison of real part of the signal')
box on;

figure;hold on;
plot(t(1:10000),imag(dataout(1:10000)),'LineWidth', 1)
plot(t(1:10000),imag(dataout_vpi(1:10000)),'LineWidth', 1)
xlabel('Time/s')
ylabel('Amplitude / V')
legend('matlab','vpi','FontName', 'Arial','FontSize',14,'Location', 'best')
set(gca, 'FontName', 'Arial', 'FontSize', 12);
% 设置X轴标签字体属性
set(get(gca, 'XLabel'), 'FontName', 'Arial', 'FontSize', 14);
% 设置Y轴标签字体属性
set(get(gca, 'YLabel'), 'FontName', 'Arial', 'FontSize', 14);
title('Comparison of imaginary part of the signal')
box on;