function [ RI ] = ComputeRI( N, r, p, k )
%Compute a approximation of an infinite impulse response of length 2*N-1.
%Only the causal number of taps is required.
% ComputeRI 函数通过对滤波器的留数、极点和直接项进行处理，
% 分别计算因果项、直接项和非因果项的脉冲响应，最终合并得到一个长度为 2*N - 1 的近似无限脉冲响应。
% 该函数在数字滤波器设计和信号处理中可用于根据部分分式展开的结果计算滤波器的脉冲响应。

[~,test]=unique(r);
if length(test) < length(r)
display('High order roots')
end

RI=zeros(1,N);
%Causal terms 因果项
index=find(abs(p)<1);
if length(index>0)
RI=sum(diag(r(index))*repmat(p(index), [1 N]).^repmat((0:N-1), [length(r(index)) 1]), 1);
end
%direct termes  直接项
if length(k)>0
RI(1:length(k))=RI(1:length(k))+k;
end
%Non Causal terms  非因果项
index=find(abs(p)>1);
if length(index)>0
RIneg=sum(diag(r(index))*(- (repmat(p(index), [1 N-1]).^repmat((-N+1:-1), [length(r(index)) 1]))),1);
else
RIneg=zeros(1,N-1);    
end
RI=[RIneg,RI];

end

