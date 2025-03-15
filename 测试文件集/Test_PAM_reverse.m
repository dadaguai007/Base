clc;clear;close all;
addpath("D:\BIT_PhD\Base_Code\Codebase_using")
% PAM transimision

SpS = 2;
Rs  = 10e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;

Vpi = 2;
Vb = -Vpi/2;


Pi_dBm = 10;
% Pi_dBm=[-10:10];
Pi = 10^(Pi_dBm/10)*1e-3; %W


param=struct();
% param.Ltotal=[0:20];
param.Ltotal = 40; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
param.Fs=Fs;


paramPD=struct();
paramPD.B =Rs;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;

%linearChannel
paramCh=struct();
paramCh.Fs = param.Fs;
paramCh.L = param.Ltotal;
paramCh.D = 16;
paramCh.Fc  = 193.1e12;
paramCh.alpha = 0.2;

paramEDFA = struct();
paramEDFA.G = paramCh.alpha*paramCh.L;
paramEDFA.NF = 4.5  ;
paramEDFA.Fc = paramCh.Fc;
paramEDFA.Fs = Fs;
paramEDFA.type='noise';

%PAM
M=4;
data_2bit=randi([0,1],log2(M),80000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;

% Mapeia bits para pulsos elétricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);

% MDRC 码
res1=[];
N=20;
for j=1:length(symbTx)/N
    x=symbTx((j-1)*N+1:j*N);
    % 计算RDS
    for i=1:N
        if i==1
            z(1)=x(1);
        else
            z(i)=z(i-1)+x(i);
        end
    end
    % 确定反转点
    z_k(j)=(z(N)/2);
    [~,k_index(j)]=min(abs(z-z_k(j)));

    % 进行反转
    res=[x(1:k_index(j)),-x(k_index(j)+1:end)];

    % 排列传输信号
    res1=[res1,res];
end


% Upsampling
symbolsUp = upsample(res1, SpS);

% %Pulso
hsqrt = rcosdesign(0.01,256,SpS,'sqrt');
hsqrt=hsqrt./max(abs(hsqrt));
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');


for i=1:length(sigTx)
    sigTx_reverse(i)=(-1).^i*sigTx(i);
end

mon_ESA(sigTx,Fs);