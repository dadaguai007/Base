% Test LUT
clc;clear;close all;

SpS = 6;
Rs  = 20e9;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;

lmbd =1550e-9;
Pi_dBm = 10;

%MZM
Vpi = 2;
Vb = -Vpi/2;
Pi = 10^(Pi_dBm/10)*1e-3; %W
rng(10)
%PAM
M=4;
data_2bit=randi([0,1],log2(M),8000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;

% Mapeia bits para pulsos eletricos
symbTx = pammod(symbols,M,0,'gray');
symbTx = pnorm(symbTx);

% Upsampling
symbolsUp = upsample(symbTx, SpS);


% %Pulso
hsqrt = rcosdesign(0.01,256,SpS,'sqrt');  
% % pulse shaping
sigTx=conv(symbolsUp,hsqrt,'same');

ref_sym_c = unique(symbTx);



% channel paramerent
hch = [0.207, 0.815, 0.207];
hch_up = upsample(hch, SpS);


% channel respon
sigCh = firFilter(hch_up, sigTx);
sigCh = pnorm(sigCh);
% sigCh is output
err = (sigTx - sigCh).';

figure;hold on;
plot(sigTx)
plot(sigCh)

close all;
mem_len=3;
 idx_data = ones(size(sigTx));
%  err(idx_data)
%  idx = find(idx_data>0)-floor(3/2);
     idx = find(idx_data>0); 
    shift=floor(mem_len/2);
N=3;
L = M^N;
%模式索引
Idx=1:L;
% 将索引生成矩阵形式
Idx = reshape(Idx, (repmat(M, 1, N)));
% 将模式数组进行转置
for i=1:(N+1)
Idx(:,:,i)=Idx(:,:,i).';
end


% 得到的矩阵是，发送端每个符号，与参考符号的差值，每个符号对应一列
 k=abs(bsxfun(@minus, sigTx(:).', ref_sym_c(:)));
 [~, sig_idx]=min(k);
 N=3;
sig_rwin = rolling_window(sig_idx, N, 'true');

sig_rwin=sig_rwin.';
% 第一行为哪一个矩阵，第二行为第几行，第三行为第几列
pattern_idx=zeros(length(sig_rwin),1);
for i=1:length(sig_rwin)
pattern_idx(i) = Idx(sig_rwin(2,i),sig_rwin(3,i),sig_rwin(1,i));
end
pattern_idx=pattern_idx(idx);
pattern_idx=circshift(pattern_idx,shift); 

ea = cal_lut_avg(err, pattern_idx, pattern_idx, L);
Error=ea(pattern_idx).';
sigTx_LUT=sigTx+Error;
plot_spectrum(sigTx_LUT,Fs);
plot_spectrum(sigTx,Fs);

sigCh_LUT = firFilter(sigTx_LUT, sigTx);
sigCh_LUT = pnorm(sigCh_LUT);
plot_spectrum(sigCh_LUT,Fs);
plot_spectrum(sigCh,Fs);
[LUT, idx_I, idx_Q] = cal_lut(sigTx, sigCh, symbTx, mem_len);
