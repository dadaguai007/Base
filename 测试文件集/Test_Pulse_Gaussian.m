clc;clear;
T0 = [0.2, 0.1, 0.05];
A0 = 1;
fs = 100;
NFFT = 1024;
t = -0.5:1/fs:0.5-1/fs;
f = fs/NFFT * (-NFFT/2:NFFT/2-1);
L = length(t);

for i = 1:length(T0)
    
    x = GaussPulse(A0,T0(i),t);
    figure(1);hold on;
    plot(t,x)
    
    Xf=fftshift(fft(x,NFFT));
    figure(2);hold on;
    plot(f,abs(Xf)/L)

    xa = FGaussPulse(A0,T0(i),f);
    figure(3);hold on;
    plot(f,xa)
end