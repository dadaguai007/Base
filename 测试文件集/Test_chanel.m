clc;clear; close all;
% 预均衡
% 包括均衡处理后的图像处理，MSE的显示
SpS = 4 ;           
Rs  = 10e6;         
Ts  = 1/Rs  ;        
Fs  = 1/(Ts/SpS);    
Ta  = 1/Fs       ;   

SNRdB = 30;
% channel paramerent
hch = [0.207, 0.815, 0.207];
hch_up = upsample(hch, SpS);

Ntaps = 31;
k = 0.005;

% rand
bits = randi([0, 1], 1, 10000);

% Mapeia bits para pulsos elétricos
symbTx = 2 * bits - 1;
symbTx = pnorm(symbTx);

% upsampling
symbolsUp = upsample(symbTx, SpS);

%pulse
% pulse = pulseShape('nrz', SpS);
% pulse = pulseShape('rrc', SpS,  2048, 0.01,Ts);
% pulse = pulse./ max(abs(pulse));
% sigTx = firFilter(pulse, symbolsUp);

hsqrt = rcosdesign(0.01,2048,SpS,'sqrt');  
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');
sigTx = pnorm(sigTx);

% channel respon
% sigCh = firFilter(hch_up, sigTx);

sigCh=filter(hch_up, 1,sigTx);

% sigCh = conv(hch_up, sigTx);
% sigCh=sigCh(1:end-length(hch_up)+1);

sigCh = pnorm(sigCh);


%  gaussiano noise
sigma2 = 1/(10.^(SNRdB/10));
noise=gaussianRealNoise([length(sigTx),1], sigma2).';
%match
sigRx1=sigCh+noise;
% sigRx1 = firFilter(pulse, sigRx1);
sigRx1=conv(sigRx1,hsqrt,'same');
sigRx = pnorm(sigRx1);
% sigRx = sigRx1;
eyediagram(sigRx,2*SpS);
title('接收的信号眼图')


% downsampling
y = downsample(sigRx,SpS);
% y = pnorm(y);

% equalizador LMS
[symbRx, h_eq, squaredError ]= LMS(y, symbTx, Ntaps, k);
symbRx = pnorm(symbRx);

eyediagram(symbRx,2*SpS)

% detector ML
% bitsRx = demodulateGray(np.sqrt(Es)*symbRx, M, constType)
if 0
sigLMS = firFilter(upsample(h_eq, SpS), sigRx1);
eyediagram(pnorm(sigLMS),2*SpS)
end

% 测设mse
L = 100;
MA_SE = conv(squaredError,ones(1, L)/L, 'same');
figure;
semilogy(MA_SE)
title('滑动窗口的MSE')

figure;
semilogy(squaredError)
title('MSE')


% 滤波器的响应
x=1:Ntaps;
figure;
stem(x,h_eq)
title('滤波器抽头响应')

%信号的分布
figure;hold on
 plot(y,'.')
 plot(symbRx,'k.')
plot(symbTx,'.')
legend('接收信号','均衡后信号','发送信号')


% 频域响应
fs = Rs; % 采样率等于符号速率
N=length(y);
f=(0:1/N:0.5-1/N);
[h, w] = freqz(h_eq, 1,f,1); % 计算数字滤波器的频率响应
%计算延时
groupDelay = grpdelay(h_eq,1, 1);
%频域延时
Hf = h.*exp(1j*2*pi*f*groupDelay);

% 将角频率 ω 转换为实际频率 x（单位：Hz）
x = w * fs / (2 * pi);
% 计算频率响应的幅度响应，并转换为分贝（dB）
Y= 20 * log10(abs(Hf));

[H_ch, w] = freqz(hch, 1,f,1); % 计算数字滤波器的频率响应
y_ch = 20 * log10(abs(H_ch));
figure;hold on;
plot(x/1e6,y_ch)
plot(x/1e6,Y)

%信道反函数与抽头响应对比
figure;hold on;
plot(abs(1./H_ch))
plot(abs(h))


% 转换为dB轴
y_ch_reverse = 20 * log10(abs(1./H_ch));
figure;hold on;
plot(x/1e6,y_ch_reverse)
plot(x/1e6,Y)




%% fifter函数 应用发射端 使用预均衡 
hch_pre = upsample(h_eq, SpS);
sigCh_revese = firFilter(hch_pre, sigTx);
sigCh_revese = pnorm(sigCh_revese);
plot_spectrum(sigCh_revese,Fs);
% 无延时
figure;
plot(sigCh_revese)

sigCh_pre = firFilter(hch_up, sigCh_revese);
sigCh_pre = pnorm(sigCh_pre);
plot_spectrum(sigCh_pre,Fs);

eyediagram(sigCh_pre,2*SpS)

%% 应用于发射端
hch_pre = upsample(h_eq, SpS);
% 使用滤波器
y_pre = filter(hch_pre, 1, sigTx);
figure;
plot(y_pre)
plot_spectrum(y_pre,Fs);

eyediagram(y_pre,2*SpS);

sigCh_pre = filter(hch_up, 1,y_pre);
sigCh_pre = pnorm(sigCh_pre);
plot_spectrum(sigCh_pre,Fs);

eyediagram(sigCh_pre,2*SpS)


%%
% 存在bug,并没有延时
% 计算并抵消群时延 ，输出延时多少采样点
delay= grpdelay(hch_pre,1,1);
%去除延时点直接删除
y(1:floor(delay))=[];
%不取
y_no_delay = y(floor(delay)+1:end);

plot_spectrum(y,Fs);


eyediagram(sigCh_pre,2*SpS)
%% 逆卷积的求取
L=512;
delta=[1 zeros(1,L-1)]; %生成长度为L的单位抽样序列
[p,q]=deconv(h,delta);
figure;
stem(q)


%% 
% 均衡器响应
L=512;
delta=[1 zeros(1,L-1)]; %生成长度为L的单位抽样序列
equalizer_response = filter(h_eq, 1, delta);
subplot(2,1,1)
plot(h_eq)
subplot(2,1,2)
plot(equalizer_response)

% 信道响应
estimated_channel_response = filter(hch, 1,delta);
figure;
plot(estimated_channel_response)
%% 逆卷积求取还是有问题
% [p,q]=deconv(estimated_channel_response,h_eq);
% figure;
% plot(q)
% [p,q]= deconv(equalizer_response, estimated_channel_response);
% figure;
% plot(q)