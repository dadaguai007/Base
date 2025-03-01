function Plot_signal_equ(rec,symbRx,symbTx)

%信号的分布
figure;hold on
plot(rec,'.')
plot(symbRx,'k.')
plot(symbTx,'.')
legend('接收信号','均衡后信号','发送信号')


end