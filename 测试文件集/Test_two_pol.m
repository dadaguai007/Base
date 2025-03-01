% Test two_pol_filter
% test_plot_filter/trsonform
clc;clear;
fs=10e9;
Rs=1e9;

verbose=1;
N=1000;
filt=two_pol_filter(Rs,fs,verbose);
% 绘制频域响应
[f, t] = freq_time_set(N, fs);
H=filt.H(f/fs);
figure;
plot(real(H))
% test_H_delay
% [H1, delay]=H_grpdelay(H, f);