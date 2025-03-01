clc;clear;close all;

% pam or qam
type='pam';
SpS = 16;
Rs  = 60e9;
Ts  = 1/Rs;
Fs  = 1/(Ts/SpS);
Ta  = 1/Fs;
rolloff = 0.001;

Pi_dBm=10;
%MZM
Vpi = 2;
Vb = -Vpi/2;
Pi = 10^(Pi_dBm/10)*1e-3; %W

%fiber
param=struct();
param.Ltotal = 40; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
% param.amp='none';
param.Fs=Fs;

%PD
paramPD=struct();
paramPD.B =Rs;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;


rng(1)

M=2;
% rand
bits = randi([0, 1], 1, 10000*log2(M));
%预编码
c_bits=zeros(size(bits,1));
for i = 1:length(bits)
    if i==1
        c_bits(i)=bits(i);
    else
        %异或运算
        c_bits(i) =xor(bits(i),c_bits(i-1));
    end
end
c_bits=reshape(c_bits,log2(M),[]);
if strcmp(type,'pam')
    symbols = 2.^(0:log2(M)-1)*c_bits;
    symbTx = pammod(symbols,M,0,'gray');
    symbTx = pnorm(symbTx);
elseif strcmp(type,'qam')
        symbTx=qammod(c_bits,M,'InputType','bit','UnitAveragePower',1) ;
%     symbols = 2.^(0:log2(M)-1)*c_bits;
%     symbTx=pskmod(symbols,M,pi/4,'gray') ;
%     symbTx = pnorm(symbTx);
end

symbolsUp = upsample(symbTx, SpS);

reverse='False';
pulseDoub = Duob(SpS, 2048, reverse, rolloff);
pulse = rcosdesign(rolloff,2048,SpS,'sqrt');
% pulse = pulseShape('rrc', SpS, 2048, rolloff);
% pulse = pulse/max(abs(pulse));
% sigTx = firFilter(pulseDoub, symbolsUp);
sigTx=conv(pulseDoub,symbolsUp,'same');
t = (0:length(symbTx)-1) * (Ta / 1e-9);

%星座图
scatterplot(downsample(sigTx,SpS))
%眼图
eyediagram(sigTx,4*SpS)

%mzm
Ai     = sqrt(Pi);
% 放大倍数
Amp=1;
sigTxo = mzm(Ai, Amp*sigTx, Vb,Vpi);


%ssfo
sigRxo=ssfm(sigTxo,param);
%pd
ipd = pd(sigRxo, paramPD);
%match
Ipd=conv(pulse,ipd,'same');
Ipd = pnorm(Ipd);
idX = 1:1024;
figure
plot(t(idX),Ipd(idX),LineWidth=1)
eyediagram(Ipd,2*SpS)
title('pd-after-fiber')
plot_spectrum(Ipd,Fs);

% 功率选择性衰落
[powershift,fshift,p1,fnew,p2,f] = FFT(Ipd,Fs);
figure;
plot(fnew,20*log(p1))
% 频谱表示
F=abs(fftshift(fft(Ipd)));
figure;
plot(f,20*log10(F))



