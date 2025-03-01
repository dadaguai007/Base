clc;clear;

%周期图法
fs = 1000;                   %采样频率1000Hz
N = 50000;                   %采样50000个点
f = (0:(N-1))*fs/N;          %功率谱的横坐标

[b,a] = cheby1(5,5,0.1);     %生成低通滤波器
in = rand(1,N);              %生成随机信号
out = filter(b,a,in);        %随机信号通过低通滤波器
window = hanning(N)';        %计算窗函数的能量（第一步）
winout = out.*window;        %计算窗函数的能量（第二步）
fout = abs(fft(winout,N)).^2;%对每一点求fft并取模方
U = sum(window.*window);     %计算窗函数的能量（第三步）
%功率谱单边的话，需要乘2。
f1out = 2*fout/U;
psd1 = 10*log10(abs(f1out)); %取对数
figure
subplot(2,1,1)
plot(f(1:25000),psd1(1:25000))
xlabel('Frequency, Hz')
ylabel('周期图法：功率谱')

subplot(2,1,2)
plot(f(1:5000),psd1(1:5000))
grid; axis([0 100 -70 10]);

xlabel('Frequency, Hz')
ylabel('周期图法：前五千个点的功率谱')



%分段周期图法

fs = 1000;                   %采样频率1000Hz
N = 50000;                   %采样50000个点

[b,a] = cheby1(5,5,0.1);     %生成低通滤波器
in = rand(1,N);              %生成随机信号
out = filter(b,a,in);        %随机信号通过低通滤波器
K = 25;                      %把50000个点分为25段
M = N/K;                     %M为每一段的点数
fK= (0:(M-1))*fs/M;          %频率横坐标
d = zeros(1,M);
psdk = zeros(1,M);
window = hanning(M)';        %计算窗能量
U = sum(window.*window);     %计算窗能量
for k = 1:K                 %嵌套for循环结构，实现分段计算fft
    for j = 1:M
        index = (k-1)*M+j;
        d(j) = out(index);
    end
    dwin = d.*window;
    psdk = (abs(fft(dwin,M)).^2)/U + psdk;
end
psd2 = 10*log10(psdk/K);
figure
plot(fK(1:250),psd2(1:250))
grid;axis([0 100 -70 10]);
xlabel('Frequency ,Hz')
ylabel('分段周期图法的功率谱：250个点')


%psd
fs = 1000;                   %采样频率1000Hz
N = 50000;                   %采样50000个点
f = (0:(N-1))*fs/N;          %PSD的横坐标
[b,a] = cheby1(5,5,0.1);     %生成低通滤波器
in = rand(1,N);              %生成随机信号
out = filter(b,a,in);        %随机信号通过低通滤波器
K = 25;                      %把50000个点分为25段
M = N/K;                     %M为每一段的点数

[pxx,fK] = pwelch(out,hanning(M),0,f,fs);
psd2 = 10*log10(pxx);
figure
plot(fK(1:5000),psd2(1:5000))
grid;
xlabel('Frequency ,Hz')
ylabel('pwelch函数估计的PSD')
