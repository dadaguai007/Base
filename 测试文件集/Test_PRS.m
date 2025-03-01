% Test 双二进制
% PR响应大致有三种形式


clc,clear,close all;
dt=0.01;
Ts=1;
t=-5*Ts:dt:5*Ts;
N=Ts/dt;
eye_num=5;%眼图数量
N_data=1000;%码元数
% 产生双极性数字信号，+1，-1 序列
d=sign(randn(1,N_data)); %randn 随机生成数字,sign 函数进行一个正负判断，正数为+1，负数为-1
% dd=sigexpand(d,N); %将输入的序列扩成间隔为 N-1 个 0 的序列
% 部分响应系统冲击响应

% The First 各自平移Ts/2

ht=sinc((t+Ts/2)./Ts);
ht1=sinc((t-Ts/2)./Ts);


gt=ht+ht1;%另一种产生 gt 的方法，即 ht 与 ht 延迟一个周期 Ts 后的叠加
% st=conv(dd,gt);
tt=-5*Ts:dt:(N_data+5)*N*dt-dt; %设置采样时间
% 比较观察 t 对 h(t)和 g(t)影响
figure;
plot(t,ht)
hold on
plot(t,gt,'r--')
hold on
grid;
legend('h(t)','g(t)');
% gtext('图 3-2 t 对 g(t)影响')


% The second term
ht=sinc((t)./Ts);
ht1=sinc((t-Ts)./Ts);


gt=ht+ht1;
% tt=-5*Ts:dt:(N_data+5)*N*dt-dt; %设置采样时间
% 比较观察 t 对 h(t)和 g(t)影响
figure;
plot(t,ht)
hold on
plot(t,gt,'r--')
hold on
grid;
legend('h(t)','g(t)');


% The third term
ht=sinc((t+Ts)./Ts);
ht1=sinc((t-Ts)./Ts);

gt=ht-ht1;
% tt=-5*Ts:dt:(N_data+5)*N*dt-dt; %设置采样时间
% 比较观察 t 对 h(t)和 g(t)影响
figure;
plot(t,ht)
hold on
plot(t,gt,'r--')
hold on
grid;
legend('h(t)','g(t)');

% 这种双二进制的频域响应以后再补


%%

%研究部分响应的特性
%部分响应通过引入有规律的码间干扰，在接收端进行消除，在频带利用率最大化的情况下实现等效的无码间干扰

%---研究第I类部分响应的时域和频域响应----
clc,clear,close all;
Tb = 1; %码元长度
fs = 103.7; %采样频率  %这里是个迷，因为和第9行的那个eps，需要调参才能使得图形完整
dt = Tb/fs; %采样时间间隔
t = -3.5*Tb:dt:3.5*Tb; %时域的时间区间
gt = 4/pi*cos(pi*t/Tb)./(1-4*t.^2/Tb^2+eps); %第I类部分响应的时域表达式

figure('NumberTitle', 'off', 'Name','第I类部分响应信号特性研究');
subplot(2,1,1);
plot(t,gt);
axis([-3.5, 3.5 -0.5 1.5]);
xlabel('时间t');
ylabel('第I类部分响应的时域信号');
grid on;

f = -2/Tb:0.001:2/Tb; %画图时显示的频率范围,为了显示效果，这里我们不采用时域信号推出的频率分辨率,自定义分辨率0.001
Gf = zeros(1,length(f));
%----第I类部分响应的频域表达式----
for j = 1:length(f)
    if abs(f(j))>1/2*Tb
        Gf(j)=0;
    else
        Gf(j)=2*Tb*cos(pi*f(j)*Tb);
    end
end
subplot(2,1,2);
plot(f,Gf); %可以看到，部分响应的频带利用率与理想低通信道相同
axis([-2, 2 -0.5 2.5]);
xlabel('频率f');
ylabel('第I类部分响应的频域特性');
grid on;

%---但是直接判决会出现“差错传播现象”，故下面进行“预编码-相关编码-模2判决”的仿真
ak = [1 0 1 1 0 0 0 1 0 1 1];
b_ini = 0; %差分码的初始码元，需要指定一个初始码元
b_k =[b_ini zeros(1,length(ak))];
%计算差分码（预编码）
for i = 1:length(ak)
    b_k(i+1) = xor(ak(i),b_k(i));
end

Ck = b_k(1:end-1) + b_k(2:end); %相关编码
rec_a = mod(Ck,2); %模2判决
disp(rec_a);
%%
clc;clear;close all;
% 创建一个时间范围
t = -10:0.01:10-0.01;

% 创建一个脉冲信号（例如，矩形波）
pulse_signal = square(t);

% 创建 sinc 函数滤波器
sinc_filter = sinc(t);
% sinc_filter = sinc(t/0.01);
% sinc_filter=rcosdesign(0,length(t),1,'normal');
% 卷积脉冲信号和 sinc 函数滤波器
filtered_signal = conv(pulse_signal, sinc_filter, 'same');

% 绘制结果
figure;
subplot(3,1,1);
plot(t, pulse_signal);
title('原始脉冲信号');

subplot(3,1,2);
plot(t,sinc_filter);
title('sinc 函数滤波器');

subplot(3,1,3);
plot(t, filtered_signal);
title('滤波后的信号');
