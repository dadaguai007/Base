% Test Besslf_filter
% test_plot_filter/trsonform
clc;clear;
fs=10e9;
Rs=5e9;
fcnorm=0.1*Rs/(fs/2);
order=5;
verbose=1;
N=1000;
filt=Besslf_filter(order,fcnorm,verbose);
% 绘制频域响应
[f, t] = freq_time_set(N, fs);
H=filt.H(f/fs);
H_test=filt.H_test(f/fs);
figure;
plot(real(H))
% test_H_delay
[H1, delay]=H_grpdelay(H_test, f);
Delay=filt.grpdelay;