close all;clc;clear;

fb=10e9;
sps=4;
fs=fb*sps;

bits = randi([0, 1], 1, 10000);
n = 0:length(bits)-1;

% Mapeia bits para pulsos elétricos
symbTx = 2 * bits - 1;
symbTx = pnorm(symbTx);

% 扩频倍数
spreading_factor = 10;

% Upsampling
symbolsUp = upsample(symbTx,sps);

%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(0.2,2^10,sps,'sqrt');

bpsk_data=conv(symbolsUp,hsqrt,'same');

figure;
mon_ESA(bpsk_data,fs);

% 将原始信号扩展为10倍
expanded_bpsk_signal = repmat(bpsk_data.', spreading_factor, 1);


% 生成PN序列
sr=[1 -1 1 -1]; % 初始化PN序列生成器
pn1=[]; % 初始化PN序列
for i=1:length(bits)
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

% Upsampling
% pnUp1 = upsample(pn1,sps);
pnUp1=resample(pn1,sps,1);

pnUp1=pnUp1.';
% 生成CDMA信号
cdmaData=expanded_bpsk_signal.*pnUp1; % 将BPSK信号与PN序列相乘

figure;
mon_ESA(cdmaData,fs*10);
% 接收端解扩
rx=cdmaData.*pnUp1;
