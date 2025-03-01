clc;clear;close all;

SpS = 6;
Rs  = 10e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;

%OOK
% rand
bits = randi([0, 1], 1, 10000);
n = 0:length(bits)-1;

% Mapeia bits para pulsos elétricos
symbTx = 2 * bits - 1;
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);

% Pulso NRZ típico
pulse = pulseShape('nrz', SpS);
pulse = pulse./ max(abs(pulse));
% filter pulse
sigTx  = firFilter(pulse, symbolsUp);
% Plota sinal
t = (0:length(symbTx)-1) * (Ta / 1e-9);
idX = 1:1024;
figure;
plot(t(idX), sigTx(idX),LineWidth=1)
ylim([-1.5,1.5])
%psd
plot_spectrum(sigTx,Fs)

% awgn
N=1:15;
for i = 1:length(N)
    Y(i,:)=awgn(sigTx,N(i),'measured');
end

% hard judgement and calc ber
for j= 1:length(N)
    hard_y = Y(j,:);
    for i= 1: length(hard_y)
        if hard_y(i)>0
            Out_Y(i)=1;
        elseif hard_y(i)<0
            Out_Y(i)=-1;
        end
    end
    Out_Y=downsample(Out_Y,SpS);
    bitSeq(j,:)=(Out_Y+1)/2;
    [ber(j),nErr(j)] = CalcBER(bitSeq(j,:),bits);
end

figure;
plot(N,ber,LineWidth=1)
xlabel('SNR/dB')


%% MI
M=2;
SNRs=[5,6,7,8,9,10];
chi=[-1,1];
chi=pnorm(chi);
X=symbTx;
SNR_power = 10.^(SNRs/10);
for idx=1:length(SNRs)
    Y = awgn_channel(X, SNRs(idx));
    % Y=awgn(X, SNRs(idx),'measured');
    pxy = Calc_px_given_y(X, Y, mean(abs(X).^2)/SNR_power(idx), chi);
    mi(idx) = log2(M)-(-mean(log2(pxy)));
end
awgn_Captacal=0.5 * log2(1+ SNR_power);
figure;
plot(SNRs,mi)
hold on
plot(SNRs,awgn_Captacal)



% berplot = BERPlot();
% if 1
%     sps=[1:1:10];
%     berplot.plot(sps,ber,1,1); hold on;
%     set(gcf, 'Position', [0, 0, 480, 400]);
%     set(gca, 'XTick', sps, 'XTickLabel', sps);
%     box on;
%     set(gca, 'LineWidth', 2);
%     set(gca, 'FontName', 'Arial', 'FontSize', 10);
% end