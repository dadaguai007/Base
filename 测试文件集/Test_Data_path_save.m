% Data_path_save

% 文件存储
% 生成路径
dir_up='D:\001-处理中\相干算法\Data';
bias_name='1';
tx_amp_idx='1';
data_path = dir_up+"/train_set/mseq/bias"+bias_name+"mA/tx_amp"+tx_amp_idx;
if(~exist(data_path,'dir'))
    mkdir(char(data_path));
end
save(data_path+"/mseq18.mat","Rx","txG","tx_amp","bias")

save(data_path+"/pam"+M+"_"+b_idx+".mat","Rx","txG","tx_amp","bias")
%
origin_rate_tmp=5e9;
M=64;
t = datetime('now');
save_path = "vol_save/"+t.Year+"."+t.Month+"."+t.Day+"/"+origin_rate_tmp/1e6+"M"+"/"+M+"pam";
