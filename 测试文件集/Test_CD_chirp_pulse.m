clc;clear;

fs = 256;                   % Taxa de amostragem
NFFT = 256;                 % Número de pontos da transformada
t = -0.5:1/fs:0.5-1/fs;      % Duração da janela do pulso
L = length(t);              % Comprimento da janela temporal
f = fs/NFFT * (-NFFT/2:NFFT/2-1);
beta2 = 20e-6;              %  GVD
z = 250;                    % Distance
T0 = 0.05;                  % Parâmetro meia largura @1/e
A0 = 1;                     % Amplitude do pulso

H = exp(1i * 0.5 * beta2 * z * (2 * pi * f).^2);

figure;
plot(f/fs, abs(H), f/fs, angle(H));


x  = GaussPulse(A0, T0, t) ;
X  = FGaussPulse(A0, T0, f) ;
XH = X.*H ;
xh  = ifftshift(ifft(XH));
%解析表达式计算的频域传输后的信号在时间域的波形：
xhA = (T0 / sqrt(T0.^2 - 1i * beta2 * z)) * exp(-(t.^2) / (2 * (T0.^2 - 1i * beta2 * z)));

%光纤的色散长度（dispersion length）
LD = (T0^2)/abs(beta2);
fprintf('CD distance %f\n',LD);

figure;hold on ;
plot(f,real(XH))
plot(f,imag(XH))
plot(f,abs(XH))
xlim([-fs/2,fs/2])

figure;hold on ;
plot(t,abs(xh*L))
plot(t,abs(xhA),'--')

figure;hold on;
plot(t,x)
plot(t,abs(xhA))

% chirp CD

C=2;
xC = CGaussPulse(A0,T0,t,C);
XC = CFGaussPulse(A0,T0,f,C);
XHC = XC.*H;
xhC=ifftshift(ifft(XHC));

Q=1+(C-1i)*beta2*z/T0.^2;
xhAC = (A0/sqrt(Q)) * exp(-0.5*((1+1i*C)/Q) * (t/T0).^2); 


figure;hold on ;
plot(f,real(XHC))
plot(f,imag(XHC))
plot(f,abs(XHC))
xlim([-fs/2,fs/2])


figure;hold on ;
plot(t,abs(xhC*L))
plot(t,abs(xhAC),'--')

figure;hold on;
plot(t,abs(xC))
plot(t,abs(xhAC))