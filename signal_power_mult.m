function power = signal_power_mult(x) 
% Calculate the total power of x. % 
%使用于 多个模式 和 多通道 的 平均功率 求取 
power = sum(mean(abs(x).^2));

end