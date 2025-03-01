clc;clear all;close all;
% 这是雷哥写的模拟滤波器，想看的时候再看吧
% m是设置的延时时间 
% generate PRBS  sequence
p = 15; %prbs 
nn = 2;
T = 1e-10 ; %10GB/s信号的速率 两个b信号间间距时间100ps=1e-10s
fs= 1/T ; %传输频率1e10hz
fn = nn*fs;%采样率20G
pl = 2^p-1; %length of prbs
num = pl; % 输出序列需要输出的长度
y = +prbs(p,num)'; %输出prbs序列


N = num; 
f = (0:(N-1))*fs/N;       %功率谱的横坐标
K = 5;                 %把num个点分为k段
M = N/K;                %M为每一段的点数

%%  查看功率谱
[pxx,fK] = pwelch(y,hanning(M),0,f,fn); %功率谱密度y1信号，window窗口长度，窗口重叠，fft数据点个数，采样率
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2) 
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
% 互相关计算
fn = nn*fs;      %采样频率
% 线性插值法上采样   
x = 1:length(y);
xi = 1:1/nn:length(y);
y1= (interp1(x,y,xi,'nearest'))';
 
[r,lags] = xcorr(y1-mean(y1),y1-mean(y1));  % 两通道相关,lag横坐标，r纵坐标
lags=lags';
[c, d] = max(abs(r));  %相关最大值c和其位数d
e =20;   % 定义所要拟合的高斯曲线最高点左右参考值
x_data = lags(d-e:d+e)-lags(d);  %选取最大值左右e个单位，并将高斯曲线移到关于0对称
y_data = r(d-e:d+e);  %选取最大值左右e个单位
 figure % prbs 互相关图
 plot(lags(d-e:d+e),r(d-e:d+e))
 hold on

 %% 切比雪夫I型滤波器
 Wp = 5.9e9/(fn/nn);
 Ws = 6.2e9/(fn/nn);
Rp=1; 
Rs=40;
[n,Wn] = cheb1ord(Wp,Ws,Rp,Rs); % Wp通带截止频率，Ws阻带截止频率，Rp通带波纹,Rs阻带衰减
[b1,a1] = cheby1(n,Rp,Wp);   % 生成低通滤波器 
[h1,w1] = freqz(b1,a1,[],fn);
out = filter(b1,a1,y);      % 随机信号通过低通滤波器
[pxx,fK] = pwelch(out,hanning(M),0,f,fn); %功率谱密度
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2)
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
title('切比雪夫I型滤波器')
% 互相关计算
fn = nn*fs;      %采样频率
% 线性插值法上采样
x = 1:length(out);
xi = 1:1/nn:length(out);
y1= (interp1(x,out,xi,'nearest'))';

[r,lags] = xcorr(y1-mean(y1),y1-mean(y1));  % 两通道相关,lag横坐标，r纵坐标
lags=lags';
[c, d] = max(abs(r));  %相关最大值c和其位数d
x_data = lags(d-e:d+e)-lags(d);  %选取最大值左右e个单位，并将高斯曲线移到关于0对称
y_data = r(d-e:d+e);  %选取最大值左右e个单位
 figure % prbs 互相关图
 plot(lags(d-e:d+e),r(d-e:d+e))
 title('切比雪夫I型滤波器')

 %%  切比雪夫II型滤波器
Wp = 5.8e9/(fn/nn);
Ws = 6.2e9/(fn/nn);
Rp=1;
Rs=40;
[n,Wn] = cheb2ord(Wp,Ws,Rp,Rs);
[b2,a2] = cheby2(n,Rs,Wn);   % 生成低通滤波器Wn截止频率
[h2,w2] = freqz(b2,a2,[],fn);
out = filter(b2,a2,y);      % 随机信号通过低通滤波器
[pxx,fK] = pwelch(out,hanning(M),0,f,fn); %功率谱密度
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2)
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
title('切比雪夫II型滤波器')
% 互相关计算
fn = nn*fs;      %采样频率
% 线性插值法上采样
x = 1:length(out);
xi = 1:1/nn:length(out);
y1= (interp1(x,out,xi,'nearest'))';

[r,lags] = xcorr(y1-mean(y1),y1-mean(y1));  % 两通道相关,lag横坐标，r纵坐标
lags=lags';
[c, d] = max(abs(r));  %相关最大值c和其位数d
x_data = lags(d-e:d+e)-lags(d);  %选取最大值左右e个单位，并将高斯曲线移到关于0对称
y_data = r(d-e:d+e);  %选取最大值左右e个单位
 figure % prbs 互相关图
 plot(lags(d-e:d+e),r(d-e:d+e))
title('切比雪夫II型滤波器')

  %% 椭圆滤波器
 Wp = 5.8e9/(fn/nn);
 Ws = 6.2e9/(fn/nn);
Rp=1;
Rs=40;
[n,Wn] = ellipord(Wp,Ws,Rp,Rs);
[be,ae] = ellip(n,Rp,Rs,Wp);   % 生成低通滤波器Wn截止频率
[he,we] = freqz(be,ae,[],fn);
out = filter(be,ae,y);      % 随机信号通过低通滤波器
[pxx,fK] = pwelch(out,hanning(M),0,f,fn); %功率谱密度
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2)
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
title('椭圆滤波器')
% 互相关计算
fn = nn*fs;      %采样频率
% 线性插值法上采样
x = 1:length(out);
xi = 1:1/nn:length(out);
y1= (interp1(x,out,xi,'nearest'))';

[r,lags] = xcorr(y1-mean(y1),y1-mean(y1));  % 两通道相关,lag横坐标，r纵坐标
lags=lags';
[c, d] = max(abs(r));  %相关最大值c和其位数d
x_data = lags(d-e:d+e)-lags(d);  %选取最大值左右e个单位，并将高斯曲线移到关于0对称
y_data = r(d-e:d+e);  %选取最大值左右e个单位
 figure % prbs 互相关图
 plot(lags(d-e:d+e),r(d-e:d+e))
title('椭圆滤波器')


  %% 巴特沃斯滤波器
Wp = 5.8e9/(fn/nn);
 Ws = 6.2e9/(fn/nn);
Rp=1;
Rs=40;
[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[bb,ab] = butter(n,Wn);   % 生成低通滤波器Wn截止频率
[hb,wb] = freqz(bb,ab,[],fn);
out = filter(bb,ab,y);      % 随机信号通过低通滤波器
[pxx,fK] = pwelch(out,hanning(M),0,f,fn); %功率谱密度
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2)
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
title('巴特沃斯滤波器')
% 互相关计算
fn = nn*fs;      %采样频率
% 线性插值法上采样
x = 1:length(out);
xi = 1:1/nn:length(out);
y1= (interp1(x,out,xi,'nearest'))';

[r,lags] = xcorr(y1-mean(y1),y1-mean(y1));  % 两通道相关,lag横坐标，r纵坐标
lags=lags';
[c, d] = max(abs(r));  %相关最大值c和其位数d
x_data = lags(d-e:d+e)-lags(d);  %选取最大值左右e个单位，并将高斯曲线移到关于0对称
y_data = r(d-e:d+e);  %选取最大值左右e个单位
 figure % prbs 互相关图
 plot(lags(d-e:d+e),r(d-e:d+e))
title('巴特沃斯滤波器')

% 添加一个高斯低通滤波器
filter_size = 21; % 滤波器尺寸
sigma = 3;        % 高斯分布的标准差

% 生成高斯低通滤波器
gaussian_filter = fspecial('gaussian', [filter_size, filter_size], sigma);

% 将高斯滤波器应用于信号
y_filtered = filter2(gaussian_filter, y, 'same');
[pxx,fK] = pwelch(y_filtered,hanning(M),0,f,fn); %功率谱密度
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2)
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
title('高斯滤波器')


 %% 理想滤波器
Wp = 5.8e9/(fn/2);
%  Ws = 6.5e9/(fs/2);
Rp=1;
Rs=40;
%[n,Wn] = buttord(Wp,Ws,Rp,Rs);
% [be,ae] = lowpass(,Wn);   % 生成低通滤波器Wn截止频率
% [he,we] = freqz(be,ae,[],fn);
out = lowpass(y,Wp);      % 随机信号通过低通滤波器
[pxx,fK] = pwelch(out,hanning(M),0,f,fn); %功率谱密度
psd2 = 10*log10(pxx); % db/hz 
figure
plot(fK,psd2)
grid;
xlabel('Frequency(Hz)')
ylabel('pwelch函数估计的PSD(dB/Hz)')
title('理想滤波器')
% 互相关计算
fn = nn*fs;      %采样频率
% 线性插值法上采样
x = 1:length(out);
xi = 1:1/nn:length(out);
y1= (interp1(x,out,xi,'nearest'))';

[r,lags] = xcorr(y1-mean(y1),y1-mean(y1));  % 两通道相关,lag横坐标，r纵坐标
lags=lags';
[c, d] = max(abs(r));  %相关最大值c和其位数d
x_data = lags(d-e:d+e)-lags(d);  %选取最大值左右e个单位，并将高斯曲线移到关于0对称
y_data = r(d-e:d+e);  %选取最大值左右e个单位
 figure % prbs 互相关图
 plot(lags(d-e:d+e),r(d-e:d+e))
title('理想滤波器')

%%滤波器响应图
figure
plot(wb,abs(hb))
hold on
plot(w1,abs(h1))
plot(w2,abs(h2))
plot(we,abs(he))
axis([0 1e10 -0.2 1.2])
grid
xlabel('Frequency(Hz)')
ylabel('Attenuation(dB)')
legend('butter','cheby1','cheby2','ellip')
figure
plot(wb,mag2db(abs(hb)))
hold on
plot(w1,mag2db(abs(h1)))
plot(w2,mag2db(abs(h2)))
plot(we,mag2db(abs(he)))
axis([0 1e10 -50 5])
grid
xlabel('Frequency(Hz)')
ylabel('Attenuation(dB)')
legend('butter','cheby1','cheby2','ellip')