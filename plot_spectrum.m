function h=plot_spectrum(signal,Fs)

win = 4096;
[Snrzf,~] = pwelch(signal,win,win/2,win,'centered','power');
f= Fs * (-0.5:1/win:0.5-1/win)/1e9;
figure
h=plot(f,10*log10(Snrzf),LineWidth=1.5);
xlabel('Hz');
ylabel('dB');

end