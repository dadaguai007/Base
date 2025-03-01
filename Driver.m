%% Driver 放大器模型
% Adjust gain to compensate for preemphasis
%首先将信号减去均值，使其均值为0，乘上放大增益，再将均值乘以调整偏置系数，使其偏置有一个新的高度
%更好的调整手段
% xt，Vgain较好赋值，Vbiasadj常用取1
function out=Driver(xt,Vgain,Vbiasadj)
out = Vgain*(xt - mean(xt)) + Vbiasadj*mean(xt);
end