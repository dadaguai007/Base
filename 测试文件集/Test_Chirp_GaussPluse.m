clc;clear;
T0 = [0.2, 0.1, 0.05];
A0 = 1;
fs = 100;
NFFT = 1024;
t = -0.5:1/fs:0.5-1/fs;
f = fs/NFFT * (-NFFT/2:NFFT/2-1);
L = length(t);
C=2;
for i = 1:length(T0)
    x = CGaussPulse(A0,T0(i),t,C);
    figure(1);hold on;
    plot(t,abs(x))

    Xf=fftshift(fft(x,NFFT));
    figure(2);hold on;
    plot(f,abs(Xf)/L)

    xa = CFGaussPulse(A0,T0(i),f,C);
    figure(3);hold on;
    plot(f,abs(xa))

    Omega=t*C/(T0(i).^2);
    figure(4);hold on;
    plot(t,Omega)
end