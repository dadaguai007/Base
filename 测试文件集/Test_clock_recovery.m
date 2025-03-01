% test clock recovery

clc;clear;close all;
M=16;
SpS = 32;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;
ppm = -200;        %Deviation of the ADC sampling rate in ppm


% signal
bits=randi([0,1],log2(M),5000);
% bit=bits(:);
% symbTx = modulateGray(bit, M, 'qam');
% symbTx=symbTx.';
symbTx=qammod(bits,M,'InputType','bit','UnitAveragePower',1) ;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% % Pulso rrc típico
pulse = pulseShape('rrc', SpS, 4001, 0.01);
pulse = pulse./ max(abs(pulse));
sigTx  = firFilter(pulse.', symbolsUp);

% %设计根升余弦脉冲成型滤波器
% hsqrt = rcosdesign(0.01,4001,SpS,'sqrt');  
% % pulse shaping
% sigTx=conv(symbolsUp,hsqrt,'same');





Fs_adc = 2*Rs*(1 + ppm/1e6);
ppm_meas = (Fs_adc-2*Rs)/(2*Rs)*1e6;
fprintf('ADC sampling rate = %.5f GS/s\n',Fs_adc/1e9);
fprintf('ADC sampling clock drift = %.2f ppm\n',ppm_meas);
paramADC = struct();
paramADC.Fs_in = Fs;
paramADC.Fs_out = Fs_adc;
paramADC.jitter_rms = 0; %400e-15
paramADC.nBits =  8;
paramADC.Vmax = max(real(sigTx));
paramADC.Vmin = min(real(sigTx));
paramADC.AAF = 'on';
paramADC.N = 1001;

sigRx = adc(sigTx, paramADC);

% figure;hold on;
% plot(real(sigRx), '-o', 'markersize', 4,'Color','r');
% plot( real(sigTx), '-o', 'markersize', 4,'Color','k');
% 
% legend('采样频率偏移','正确采样')

%clock recovery with Gardner's algorithm
paramCLKREC = struct();
paramCLKREC.isNyquist = 'True';
paramCLKREC.returnTiming = 'True';
% 需要调参
paramCLKREC.ki = 1e-6;
paramCLKREC.kp = 6e-4;
paramCLKREC.maxPPM = ppm;
[outCLK, ted_values]=Test_gardnerClockRecovery(sigRx, paramCLKREC);

scatterplot(downsample(sigRx,2))
scatterplot(downsample(outCLK,2))










% % resample signal to non-integer samples/symbol rate
% downSample = 7.99955;
% 
% derta_Fs = (Fs/downSample-Fs/8)/(Fs/8)*1e6;
% 
% fprintf('sampling clock deviation %.2f ppm\n',derta_Fs)
% 
% sigRxRef = clockSamplingInterp(sigTx, Fs, Fs/8, 0);
% 
% % # ADC input parameters
% paramADC = struct();
% paramADC.Fs_in = Fs;
% paramADC.Fs_out = Fs/downSample;
% paramADC.jitter_rms = 400e-15;
% paramADC.nBits =  8;
% paramADC.Vmax = 0.5;
% paramADC.Vmin = -0.5;
% paramADC.AAF = 'off';
% paramADC.N = 1001;
% 
% sigRx = adc(sigTx, paramADC);
% 
% %clock recovery with Gardner's algorithm
% paramCLKREC = struct();
% paramCLKREC.isNyquist = 'True';
% paramCLKREC.returnTiming = 'True';
% paramCLKREC.ki = 1e-6;
% paramCLKREC.kp = 2e-4;
% 
% [outCLK, ted_values]= gardnerClockRecovery(sigRx, paramCLKREC);
% 
% scatterplot(sigRx)
% scatterplot(outCLK)

