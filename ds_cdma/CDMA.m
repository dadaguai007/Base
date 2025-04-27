% 清除所有图形窗口和命令窗口的内容
close all;
clc;
clear;
% 定义二进制数据序列
% N=6; 
% data=input('Enter Binary Data In Form of 0 and 1 in [ ] : ');
data=[1 0 0 1 1 0 1 1];

% 创建一个名为“Message BPSK Modulation”的图形窗口
figure('Name','Message BPSK Modulation','NumberTitle','off')

% 绘制原始消息信号
subplot(2,2,1);
plot(rectpulse(data,100)); % 使用rectpulse函数将数据脉冲整形为矩形波形
axis([0 length(rectpulse(data,100)) -0.2 1.2]); % 设置坐标轴范围
title('Message Signal'); % 设置标题
xlabel('n'); % 设置x轴标签
ylabel('x(n)'); % 设置y轴标签
grid on; % 开启网格

% 将数据中的0替换为-1，以便后续调制处理
data(data(:)==0)=-1;
length_data=length(data); % 获取数据长度

% 定义载波频率和其他参数
fc1=10; 
eb=2; 
tb=1; 
T=1;

% 生成矩形脉冲调制后的消息信号
msg=rectpulse(data,100);
subplot(2,2,2);
plot(msg); % 绘制消息信号
title('Message Signal in NRZ form'); % 设置标题
xlabel('n'); % 设置x轴标签
ylabel('x(n)'); % 设置y轴标签
axis([0 100*length_data -1.2 1.2]); % 设置坐标轴范围
grid on; % 开启网格

% 生成BPSK调制信号
N=length_data;
Tb = 0.0001; % 定义码元间隔
nb=100; % 每个码元的采样点数
br = 1/Tb; % 比特率
Fc = br*10; % 载波频率
t2 = Tb/nb:Tb/nb:Tb; % 时间向量
t1=0.01:0.01:length_data; % 时间向量
bpskmod=sqrt(2/T)*cos(2*pi*fc1*t1); % 生成BPSK调制信号
bpsk_data=msg.*bpskmod; % 将消息信号与BPSK调制信号相乘

% 绘制BPSK调制信号
subplot(2,2,3)
plot(bpsk_data)
title(' BPSK signal');
xlabel('Time Period(t)');
ylabel('x(t)');
axis([0 100*length_data -2 2]);
grid on;

% 绘制BPSK信号的FFT频谱
subplot(2,2,4);
plot(real(fft(bpsk_data)));
title('FFT of BPSK signal');
xlabel('Frequency');
ylabel('PSD');
grid on;

% 生成PN序列
sr=[1 -1 1 -1]; % 初始化PN序列生成器
pn1=[]; % 初始化PN序列
for i=1:length_data
    for j=1:10 
        pn1=[pn1 sr(4)]; % 生成PN序列
        if sr (4)==sr(3) 
            temp=-1;
        else 
            temp=1;
        end
        sr(4)=sr(3); % 移位寄存器操作
        sr(3)=sr(2);
        sr(2)=sr(1);
        sr(1)=temp;
    end
end

% 创建一个名为“PN Generation and CDMA”的图形窗口
figure('Name','PN Generation and CDMA','NumberTitle','off');

% 绘制PN序列
subplot(2,2,1);
stem(pn1);
axis([0,length(pn1),-1.2,1.2])
title('PN sequence for data')
xlabel('n');
ylabel('x(n)');
grid on;

% 对PN序列进行上采样
pnupsampled1=[];
len_pn1=length(pn1);
for i=1:len_pn1
    for j=0.1:0.1:tb
        pnupsampled1=[pnupsampled1 pn1(i)];
    end
end

% 绘制上采样后的PN序列
subplot(2,2,2)
stem(pnupsampled1);
axis([0,length(pnupsampled1),-1.2,1.2])
title('PN sequence for data upsampled');
xlabel('n');
ylabel('x(n)');
grid on;

% 生成CDMA信号
subplot(2,2,3);
sigtx1=bpsk_data.*pnupsampled1; % 将BPSK信号与PN序列相乘
plot(sigtx1);
title('CDMA Signal');
xlabel('Time Period(t)');
ylabel('x(t)');

% 绘制CDMA信号的FFT频谱
subplot(2,2,4);
plot(real(fft(sigtx1)));
title('FFT of spreaded CDMA Signal');
xlabel('Frequency');
ylabel('PSD');
grid on;

% 添加高斯白噪声
sigtonoise=20; % 定义信噪比
composite_signal=awgn(sigtx1,sigtonoise);

% 创建一个名为“CDMA Reciever”的图形窗口
figure('Name','CDMA Reciever','NumberTitle','off')
subplot(2,2,1);
plot(sigtx1);
title(' Tx Signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

% 绘制加噪后的信号
subplot(2,2,2);
plot(composite_signal);
title(sprintf('Tx signal + noise\n SNR=%ddb',sigtonoise));
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

% 接收端解扩
rx=composite_signal.*pnupsampled1;
subplot(2,2,3);
plot(rx);
title('CDMA Demodulated signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

% BPSK解调
y=[];
bpskdemod=rx.*bpskmod;
for i=1:100:size(bpskdemod,2)
    y=[y trapz(t1(i:i+99),bpskdemod(i:i+99))]; % 积分解调
end
y(y(:)<=0)=-1; % 判决逻辑
y(y(:)>=0)=1;
rxdata=y;

% 绘制解调后的信号
subplot(2,2,4);
plot(rectpulse(rxdata,100)); 
axis([0 length(rectpulse(rxdata,100)) -1.2 1.2]);
title('Recieved Message Signal in NRZ');
xlabel('n');
ylabel('x(n)');
grid on;

% 将解调后的信号恢复为原始二进制数据
rxdata(rxdata(:)==-1)=0;
rxdata(rxdata(:)==1)=1;
rxmsg=rxdata;

% 创建一个名为“Diffrent SNR”的图形窗口
figure('Name','Diffrent SNR','NumberTitle','off')
subplot(3,1,1)
plot(rectpulse(rxmsg,100)); 
axis([0 length(rectpulse(rxmsg,100)) -0.2 1.2]);
title('Recieved Message Signal');
xlabel('n');
ylabel('x(n)');
grid on;

% 生成不同信噪比下的加噪信号
sigtonoise1=5;
composite_signal1=awgn(sigtx1,sigtonoise1);
subplot(3,1,2);
plot(composite_signal);
title(sprintf('Tx signal + noise\n SNR=%ddb',sigtonoise1));
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

sigtonoise2=0;
composite_signal2=awgn(sigtx1,sigtonoise2);
subplot(3,1,3);
plot(composite_signal2);
title(sprintf('Tx signal + noise\n SNR=%ddb',sigtonoise2));
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

% 绘制星座图
scatterplot(composite_signal); grid minor;
title('Constellation Diagram of BPSK with Noise')
grid on
scatterplot(bpsk_data); grid minor;
title('Constellation Diagram of BPSK')
grid on