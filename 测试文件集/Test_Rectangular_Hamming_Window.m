% 矩形窗和汉明窗的演示


x=linspace(0,100,10000);
h=zeros(10000,1);
h(1001:8000)=1;
subplot(2,2,1);
plot(h);
title('Time domain waveform(Rectangular)');
xlabel('Samples')
ylabel('Amplitude');
hold on;
%axis([0 10000 -1 2]);
w=h(1001:1060);
W=fft(w,1024);
W2=W/W(1);
W3=20*log10(abs(W2));
W4=2*[0:511]/1024;
W5=[0:1023]/1024;
subplot(2,2,2);
%plot(W4,W3(1:512));
plot(W5,W3);
title('Window characteristics(Retangular)');
xlabel('Normalized frequency')
ylabel('Amplitude/dB');
hold on;


%haming 
ham=hamming(60);
subplot(2,2,3);
plot(ham);
title('Time domain waveform(Hamming)');
xlabel('Samples')
ylabel('Amplitude');
hold on;
H=fft(ham,1024);
H1=H/H(1);
H2=20*log10(abs(H1));
H3=2*[0:511]/1024;
H5=[0:1023]/1024;
subplot(2,2,4);
%plot(W4,W3(1:512));
plot(H5,H2);
title('Window characteristics(Hamming)');
xlabel('Normalized frequency')
ylabel('Amplitude/dB');