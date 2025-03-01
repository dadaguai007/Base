clc;clear;
%16QAM
M=16;
SpS = 2;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W

%IQ
paramIQ=struct();
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;



%PD
paramPD=struct();
paramPD.B =Rs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=Fs;


%LO  
Plo_dBm  = 10;   
f_lo = 150e6 ;  
phi_lo  = 0 ;    
lw = 500e3;
Plo =10.^(Plo_dBm/10)*1e-3 ;


% signal
bits=randi([0,1],log2(M),2^13);
symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso rrc típico
pulse = pulseShape('rrc', SpS, 4096, 0.01, Ts);
pulse = pulse./ max(abs(pulse));

%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(0.1,40,SpS,'sqrt');  
% pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');
% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);




%% optical modulator
Ai  = sqrt(Pi);
Amp=0.5;
sigTxo = iqm(Ai, Amp*sigTx, paramIQ);



% 引入相位误差（例如，由于载波频率偏移）
phase_error = 2*pi*0.1*(0:length(sigTxo)-1)/length(sigTxo); % 10%的频率偏移
E_with_error = sigTxo .* exp(1i*phase_error);

% 调用相位恢复函数
[E_rec, phase] = phase_partition_16qam(E_with_error.', SpS);

% 比较恢复后的信号与原始信号
error = E - E_rec;
mse = mean(abs(error).^2); % 均方误差

% 绘制结果
figure;
subplot(2,1,1);
plot(real(E), imag(E), 'b.', real(E_rec), imag(E_rec), 'r+');
title('原始信号与恢复信号');
xlabel('实部');
ylabel('虚部');
legend('原始信号', '恢复信号');

subplot(2,1,2);
plot(phase_error, 'b', phase, 'r');
title('引入的相位误差与估计的相位');
xlabel('样本点');
ylabel('相位（弧度）');
legend('引入的相位误差', '估计的相位');

fprintf('均方误差（MSE）: %f\n', mse);
