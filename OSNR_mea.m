% osnr
% snr和OSNR转换函数
function osnr=OSNR_mea(lmbd,snr,Rs)

c=299792458;
% lmbd =1550e-9;
f=0.1e-9;
% snr=[1:10];
% Rs=10e9;

% 参考带宽
Bref=c*f/(lmbd.^2);

% osnr
osnr=Rs*snr/(2*Bref);

end