% Test_Gaussion_filter
clc;clear;
order=2;
fs=20e9;
f3dB=5e9;
fcnorm = f3dB/(fs/2);
nlength = 1000; 
verbose=1;
N=1000;
filt=Gaussian_filter(order,fcnorm,nlength,verbose);
[f, t] = freq_time_set(N, fs);
H=filt.H(f/fs);
figure;
plot(real(H))
