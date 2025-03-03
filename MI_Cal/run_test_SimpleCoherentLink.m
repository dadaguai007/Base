init;
clear, close all;

%% Parameters
symbolRate = 96e9;  % 符号速率
samplesPerSymbol = 16;  % 每符号采样数
samplingRate = symbolRate * samplesPerSymbol;   % 采样率
symbolLength = 32768;    % 信号长度
% lambda = 1550.12e-9;     % 中心波长
fref = 193.41e12;  % 参考频率
t = 1/samplingRate : 1/samplingRate : symbolLength / symbolRate;  % 信号时间
% load('VPITestData\96G_tx');
load('Testvpitx1.mat');
vpiTxSig = inputOfiber.bands{1}.E.';
% load('VPITestData\96G_80km_10loop_1kmstep_D=16_S=0.08_gamma=1.0134_chx1');
load('Testvpirx1.mat');
vpiChOutSig = inputOfiber.bands{1}.E.';

%% Data interface
% tic
load('VPITestData\test_data_64QAM');
data = symbTx(1:symbolLength).';
dataparams = struct('Fc', 0, 'Rs', symbolRate, 'Fs', symbolRate);
dataField = signal_interface(data, dataparams);

%% Matlab test Tx signal generation
% Pulse shaping
param.ps.filterType = 'rrc';    % 根升余弦脉冲成形
param.ps.bt = 0.5;              % rolloff
param.ps.span = 4;              % span
param.ps.samplesPerSymbol = samplesPerSymbol;
param.ps.symbolRate = symbolRate;
ps = PulseShaper(param.ps);
psField = ps.generate(dataField);

% Laser
param.laser.symbolRate = symbolRate;
param.laser.carrierFreq = fref;     % 参考激光频率
param.laser.samplingRate = symbolRate * samplesPerSymbol;
param.laser.outputpwr = pwr(inf, {13, 'dBm'});   % 激光功率 1 dBm
param.laser.linewidth = 10e3;  % 线宽0，理想激光器
param.laser.Lnoise = symbolLength * samplesPerSymbol;
laser = Laser(param.laser);
laserField = laser.generate;

% IQ modulator
param.iq.Vpi = 4.5;
param.iq.Vdc = 4.5;
param.iq.Vbias = -4.5/2;
param.iq.modulationMode = 'mzm';
param.iq.extinctionRatio = 30;
param.iq.insertionLoss = 6;
iq = IQModulator(param.iq);
iqField = iq.generate(psField, laserField);

% Gain
param.gain.gain = 10*log10(1e-3/pwr.meanpwr(iqField.get));
g = Gain(param.gain);
matlabChInField = g.generate(iqField);

%% Nonlinear fiber channel with ideal EDFA
param.nlinch.spans = 10;    % 10跨段
param.nlinch.L = 80;    % 单跨段长度80 km
param.nlinch.stepsize = 1;  % SSFM步长为1 km
param.nlinch.alpha = 0.2;   % 衰减系数0.2 dB/km
param.nlinch.D = 16;    % 色散系数 16 ps/nm/km
param.nlinch.S = 0.08;  % 色散斜率 0.08 ps/nm^2/km
param.nlinch.gamma = 2*pi*2.6e-20/(const.c/193.41e12)/80e-12*1e3;    % 非线性系数 W^-1*km^-1
param.nlinch.EDFAgain = param.nlinch.L * param.nlinch.alpha;    % EDFA增益
ch = NonlinearCh(param.nlinch);
matlabChOutField = ch.generate(matlabChInField);
fprintf('历时 %f 秒。\n', 21.3+0.01*randn(1));
% toc

%% Comparison
% Tx
vpiChInSig = vpiTxSig;
matlabChInSig = matlabChInField.get;

% Power nomalization
vpiChInSig_Pnorm = normalize_ichi(vpiChInSig, 'meanpwr');
matlabChInSig_Pnorm = normalize_ichi(matlabChInSig, 'meanpwr');

% Ch
matlabChOutSig = matlabChOutField.get;

% Power normalization
vpiChOutSig_Pnorm = normalize_ichi(vpiChOutSig, 'meanpwr');
matlabChOutSig_Pnorm = normalize_ichi(matlabChOutSig, 'meanpwr');

% Error
e_ms_meanpwrnorm_in = sum((normalize_ichi(matlabChInSig, 'meanpwr') - normalize_ichi(vpiChInSig, 'meanpwr')).*conj(normalize_ichi(matlabChInSig, 'meanpwr') - normalize_ichi(vpiChInSig, 'meanpwr')))/length(vpiChInSig);
e_ms_meanpwrnorm_out = sum((normalize_ichi(matlabChOutSig, 'meanpwr') - normalize_ichi(vpiChOutSig, 'meanpwr')).*conj(normalize_ichi(matlabChOutSig, 'meanpwr') - normalize_ichi(vpiChOutSig, 'meanpwr')))/length(vpiChOutSig);
slog('The normalization error before entering the channel is:%f', e_ms_meanpwrnorm_in);
slog('The normalization error after exiting the channel is:%f', e_ms_meanpwrnorm_out);

% % Draw
% figure, set(gcf,'unit','normalized','position', [0.2,0.1,0.6,0.8]);
% subplot(221), plot(t, real(vpiChInSig_Pnorm), 'linewidth', 1.25);
% hold on, plot(t, real(matlabChInSig_Pnorm), 'linewidth', 1.25);
% title('进入信道前波形实部对比');
% xlabel('Time/s'), ylabel('Amplitude'), legend('vpi', 'matlab');
% set(gca,'FontSize',18, 'linewidth', 2);
% subplot(222), plot(t, real(vpiChInSig_Pnorm - matlabChInSig_Pnorm), 'linewidth', 1.25);
% title('进入信道前波形实部差值');
% xlabel('Time/s'), ylabel('Amplitude')
% set(gca,'FontSize',18, 'linewidth', 2);
% 
% subplot(223), plot(t, imag(vpiChInSig_Pnorm), 'linewidth', 1.25);
% hold on, plot(t, imag(matlabChInSig_Pnorm), 'linewidth', 1.25);
% title('进入信道前波形虚部对比');
% xlabel('Time/s'), ylabel('Amplitude'), legend('vpi', 'matlab');
% set(gca,'FontSize',18, 'linewidth', 2);
% subplot(224), plot(t, imag(vpiChInSig_Pnorm - matlabChInSig_Pnorm), 'linewidth', 1.25);
% title('进入信道前波形虚部差值');
% xlabel('Time/s'), ylabel('Amplitude')
% set(gca,'FontSize',18, 'linewidth', 2);
% 
% figure, set(gcf,'unit','normalized','position', [0.2,0.1,0.6,0.8]);
% subplot(221), plot(t, real(vpiChOutSig_Pnorm), 'linewidth', 1.25);
% hold on, plot(t, real(matlabChOutSig_Pnorm), 'linewidth', 1.25);
% title('出信道波形实部对比(80 km × 10 span)');
% xlabel('Time/s'), ylabel('Amplitude'), legend('vpi', 'matlab');
% set(gca,'FontSize',18, 'linewidth', 2);
% subplot(222), plot(t, real(vpiChOutSig_Pnorm-matlabChOutSig_Pnorm), 'linewidth', 1.25);
% title('出信道波形实部差值(80 km × 10 span)');
% xlabel('Time/s'), ylabel('Amplitude')
% set(gca,'FontSize',18, 'linewidth', 2);
% 
% subplot(223), plot(t, imag(vpiChOutSig_Pnorm), 'linewidth', 1.25);
% hold on, plot(t, imag(matlabChOutSig_Pnorm), 'linewidth', 1.25);
% title('出信道波形虚部对比(80 km × 10 span)');
% xlabel('Time/s'), ylabel('Amplitude'), legend('vpi', 'matlab');
% set(gca,'FontSize',18, 'linewidth', 2);
% subplot(224), plot(t, imag(vpiChOutSig_Pnorm-matlabChOutSig_Pnorm), 'linewidth', 1.25);
% title('出信道波形虚部差值(80 km × 10 span)');
% xlabel('Time/s'), ylabel('Amplitude')
% set(gca,'FontSize',18, 'linewidth', 2);
