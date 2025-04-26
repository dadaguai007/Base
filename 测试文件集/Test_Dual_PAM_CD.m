% 复现论文：Optical multi-path interference mitigation for PAM4-IMDD systems using balanced coding
% 论文：Direct detection transmission of a PAM signal with power fading mitigation based on Alamouti coding and dual-drive MZM

clc;clear;close all;
addpath("D:\BIT_PhD\Base_Code\Codebase_using")
% PAM transimision

SpS = 2;
Rs  = 40e9;
Ts  = 1/Rs ;
Fs  = SpS*Rs;
Ta  = 1/Fs;

Vpi = 2;
Vb = -Vpi/2;



param=struct();
param.Ltotal = 10; %km
param.Lspan =5;
param.hz= 0.1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=Fs;





%PAM
M=4;
data_2bit=randi([0,1],log2(M),80000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;

% Mapeia bits para pulsos elétricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);

% MDRC 编码
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
symbolsUp_reverse = upsample(res1, SpS);

% %Pulso
hsqrt = rcosdesign(0.5,256,SpS,'sqrt');
hsqrt=hsqrt./max(abs(hsqrt));
% pulse shaping
sigTx=conv(symbolsUp_reverse,hsqrt,'same');


for i=1:length(sigTx)
    sigTx_reverse(i)=(-1).^i*sigTx(i);
end

% phase noise 
data_num=1e6;
lw      = 1e3;    % laser linewidth
phi_pn_lo = phaseNoise(lw, data_num, Ta);
sigLO = exp(1j * phi_pn_lo);

% Alamouti Code
symbolsUp = upsample(symbTx, SpS);
sigTx_alamouti=conv(symbolsUp,hsqrt,'same');
codedData = alamouti_code_pam(sigTx_alamouti);

% 符号相反的Alamouti Code
codedData2 = alamouti_code_pam_reverse(sigTx_alamouti);
% 测试两个向量组进行相加，
Code=[codedData,codedData2];

% MDRC码进行Alamouti编码
codedData_MDRC=alamouti_code_pam(sigTx);

% 奇数为负，偶数为正
codedData_Upper_lower = Upper_Lower_code_pam(sigLO);

% Alamouti——相加编码 
codedData3 = alamouti_code_pam_linear(sigTx_alamouti);
%LO
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi);
Pin=Ai;


% EA放大
amp_factor=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       mzm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sig_Tx=amp_factor*codedData3(1,:);
sigTxo = mzm(Pin, Sig_Tx, Vb,Vpi);

E_x=ssfm(sigTxo,param);
I_x=(E_x).*conj(E_x);
% figure;
% mon_ESA(I_x,Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     DD-mzm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vbias1=Vpi/4;
Vbias2=-Vpi/4;

VRF1=amp_factor*codedData3(1,:);
VRF2=amp_factor*codedData3(2,:);
[Eout] = MZM_DualDrive(Pin,VRF1,VRF2,Vbias1,Vbias2,Vpi,1);

Eo=ssfm(Eout,param);
I_pd= (Eo).*conj(Eo);
% I=(Eo).*conj(Eo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DD-MZM 奇偶相加                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(I_pd)/2
    I_odd(i)=I_pd(2*i-1);
    I_even(i)=I_pd(2*i);
    I_add(i)=I_odd(i)+1j*I_even(i);

end


% I_add_one=I_add(1:length(I_add)/2);
% I_add_two=I_add(length(I_add)/2+1:end);
% 
% test=I_add_one+I_add_two;
figure;
mon_ESA(I_odd,Fs);
figure;
mon_ESA(I_even,Fs);
figure;
mon_ESA(I_add,Fs);
figure;
mon_ESA(I_pd,Fs);