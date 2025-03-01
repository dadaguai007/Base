clear;clc;
T0 = 0.1;
m  = [1, 2, 10, 100];
A0 = 1;
fs = 100;
NFFT = 1024;
t = -0.5:1/fs:0.5-1/fs;
f = fs/NFFT * (-NFFT/2:NFFT/2-1);
L = length(t);
for i = 1:length(m)
    x=SuperGaussPulse(A0,T0,t,m(i));
    figure(1);hold on;
    plot(t,x)

    Xf = fftshift(fft(x,NFFT));
    figure(2);hold on;
    plot(f,abs(Xf)/L)

end