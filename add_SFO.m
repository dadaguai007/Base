function y_SFO=add_SFO(y,Nsym,fs_tx,delta_fs)
fs_rx = fs_tx + delta_fs; % 接收端采样频率（Hz）
T = Nsym/fs_tx; % 信号持续时间（s）
t_tx = 0:1/fs_tx:T-1/fs_tx; % 发射端时间向量
t_rx = 0:1/fs_rx:T-1/fs_rx; % 接收端时间向量

%插值
y_SFO = interp1(t_tx,y,t_rx,"spline"); % 对接收信号进行重采样
end