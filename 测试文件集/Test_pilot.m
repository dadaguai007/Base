
% 生成导频信号
% filter_10m为信道响应
pilot_ini = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
pilot_bpsk_ini = pilot_ini*2-1;
pilot_bpsk_tmp = conv(pilot_bpsk_ini,filter_10m);
pilot_bpsk_tmp = pilot_bpsk_tmp((length(filter_10m)+1)/2 : length(pilot_bpsk_tmp)-(length(filter_10m)-1)/2);
pilot_bpsk_tmp(pilot_bpsk_tmp <= 0) = -1;
pilot_bpsk_tmp(pilot_bpsk_tmp > 0) = 1;
pilot_bpsk = pilot_bpsk_tmp;





% 发送信号
signal_ori = [pilot_bpsk,zeros(1,zero_length),data_mpam];
pilot_bpsk_forsyn = [pilot_bpsk,zeros(1,zero_length_forsyn)];
% 上采样
pilot_upsample_forsyn = resample(pilot_bpsk_forsyn,filter_transmit,upf_transmit,dof_transmit);
signal_upsample = resample(signal_ori,filter_transmit,upf_transmit,dof_transmit);

pilot_send_forsyn = pilot_upsample_forsyn./norm(pilot_upsample_forsyn,2)*sqrt(length(pilot_upsample_forsyn))*100*1.1^(amp_inf);
signal_send_tmp = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp);
% signal_send_tmp_inf = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp_inf);

% signal_rand = rand([1,512]);
% signal_send = [signal_rand pilot_send_forsyn signal_send_tmp];
% signal_send_inf = [signal_rand pilot_send_forsyn signal_send_tmp_inf];











function[pilot] = ruo_pilot_gen(pri_poly)
n = length(pri_poly); % 多项式的长度
N = 2^n-1; % 长度
register = [zeros(1,n-1) 1];
pilot = zeros(1,N);
%通过将多项式与当前寄存器状态的元素相乘、求和并取模 2 来实现的。
% 更新寄存器状态，将寄存器的最低位输出作为 m-序列的当前位。
% 最后一位为寄存器的最低位
for i = 1:N
    newregister = mod(sum(pri_poly.*register),2);
    % 寄存器后移
    register(2:n) = register(1:(n-1));
    register(1) = newregister;
    pilot(i) = register(n);
end
end