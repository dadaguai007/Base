
% 例如生成一个AR(1)过程的有色噪声
ar_order = 1;
a = [1 -0.8]; % AR系数
n = 1000;     % 样本数量
e = randn(n, 1); % 白噪声输入
x = filter(a, 1, e); % 生成有色噪声信号

Cxx = cov(x);
L = chol(Cxx, 'lower'); % 计算下三角Cholesky因子
W = inv(L'); % Cholesky分解得到的逆矩阵作为白化滤波器
y = W.* x; % 应用白化滤波器

[autocorr_y,lag] = xcorr(y, 'unbiased');
fs=10e9;
psd_y = pwelch(y, [], [], [], fs); % 假设fs为采样频率
figure;
subplot(2,1,1); plot(lag, autocorr_y); title('Autocorrelation of Whitened Noise');
subplot(2,1,2); plot(db(psd_y)); title('Power Spectral Density of Whitened Noise');


%%
clc
close all
clearvars
w=-2*pi:(1e-2)*pi:2*pi;
z=exp(1j*w);
S=(1.04+0.2*(z+z.^(-1)))./(1.25+0.5*(z+z.^(-1)));
figure()
plot(w./pi,abs(S),'g-','LineWidth',2)
title('有色噪声功率谱密度')
xlabel('\omega \times\pi(rad)')
ylabel('S_n(\omega)')
grid on

H=(z+0.5)./(z+0.2);
So=(H.^2).*S;
figure()
plot(w./pi,abs(So),'b-','LineWidth',2)
title('白化滤波器输出噪声功率谱密度')
xlabel('\omega \times\pi(rad)')
ylabel('S_o(\omega)')
grid on

R=ifft(So);
figure()
plot([0:length(R)-1],fftshift(abs(R)),'Color',[0.0 0.8 0.9])
title('白化滤波器输出噪声的自相关函数')
xlabel('\tau')
ylabel('R(\tau)')
grid on


%%
clc;clear;close all;
% AR_Whitening
% 参数设定：
% 确定滤波器的阶数N，阶数应大于1
N=10;

% 设置绘制功率谱参数，参数详情可参见pwelch函数说明
Fs=10000; % 原始数据采样率
nwindow=100000; % 窗长
coe=1; % 数据乘系数，默认为1
nfft=100000;
overlap=0.5;
% 例如生成一个AR(1)过程的有色噪声
ar_order = 1;
a = [1 -0.8]; % AR系数
n = 1e6;     % 样本数量
e = randn(n, 1); % 白噪声输入
X = filter(a, 1, e).'; % 生成有色噪声信号
% 白化滤波
OriData=X(1,:);

[LayerNum,DataLength]=size(OriData);

FullResultMat=zeros(LayerNum,DataLength);

for i=1:LayerNum

    ChosenData=OriData(i,:);
    % 计算互相关函数 r(x)
    [Corr,lg]=xcorr(ChosenData);
    Corr(lg<0) = [];
    % CorrDisplay=Corr(1,1:N+1);

    % 由Levinsion-Durbin算法求得反射系数 K
    [A,E,K]=levinson(Corr,N);

    % 由格型结构计算白化结果
    [WhiteningResult,BackwardError]=latcfilt(K,ChosenData);

    FullResultMat(i,:)=WhiteningResult;

end



% 绘制功率谱

% 设置功率谱绘制数据区间
nbegin=1;
nend=DataLength;

x=(OriData(:,nbegin:nend)*coe)';
y=(FullResultMat(:,nbegin:nend)*coe)';
window=hanning(nwindow);

noverlap=overlap*nfft;
[Pxx,f1]=pwelch(x,window,noverlap,nfft,Fs);
[Pyy,f2]=pwelch(y,window,noverlap,nfft,Fs);
[Pee,f3]=pwelch(e,window,noverlap,nfft,Fs);
ee=sqrt(Pee*2);
xx=sqrt(Pxx*2);
yy=sqrt(Pyy*2);
figure(5)
loglog(f3,ee);
title('白噪声功率谱')
figure(1); % 图1：原数据功率谱
loglog(f1,xx);
title('原数据功率谱')
figure(4)
plot(X)
title('原数据时域图')
figure(2); % 图2：白化滤波数据功率谱
loglog(f2,yy);
title('白化滤波数据功率谱')
figure(3); plot(FullResultMat); % 图3：白化滤波时域结果
title('白化滤波时域图')