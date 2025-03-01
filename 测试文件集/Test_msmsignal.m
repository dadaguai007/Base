clear,
rng(12);
% 多符号调制
% Parameter Definition
sps = 2; % Samples per symbol
leng = 2e6 + 20*sps; % Data length
bt = 1/8; % Multi-Symbol Modulation(MSM) bandwidth limitation
span = 2; % ISI Symbol num or MSM order

% Data Interface
data = randi([0,1],leng,1);

% qpsk_modulation
xi = data(1:2:end) * 2 - 1;
xq = data(2:2:end) * 2 - 1;
qpsk_signal = xi + 1i*xq;

% upsampling
upsample_signal = upsample(qpsk_signal,sps); % Upsample

% Multi-Symbol Modulation
mod_ch = gaussdesign(bt,span,sps);
msm_signal = filter(mod_ch,1,upsample_signal);

msm_signal = msm_signal(10*sps*sps+1:end);
scatterplot(msm_signal)


