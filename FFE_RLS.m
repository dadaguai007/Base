function [ffe,en,w] = FFE_RLS(obj, xn, dn)
dn = dn(:);
%信号长度
L1 = length(xn);
L2 = length(dn);
%滤波器输出长度
n = round(L1/obj.sps);

% 抽头长度
w = zeros(obj.k1, 1);
x = zeros(obj.k1, 1);
yn = zeros(n,1);
en = zeros(n,1);
%参考点之间的拼接，输入向量从ref到结尾，参考向量从倒数参考点到结尾
x3 = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));
y_d = zeros(size(yn));

SD = eye(length(w));
for i = 1:n - 1
    x = cat(1,x(obj.sps+1:end),x3(obj.sps*i-obj.sps+1:1:obj.sps*i));
    yn(i) = x.'*w ;
    % 判决
    if i > L2
        if yn(i) > 2
            y_d(i) = 3;
        elseif yn(i) > 0
            y_d(i) = 1;
        elseif yn(i) > -2
            y_d(i) = -1;
        else
            y_d(i) = -3;
        end
        en(i) = y_d(i)- yn(i);
    else
        en(i) = dn(i)- yn(i);
    end
    %使用rls更新抽头
    SD = ( SD - ((SD * x) * (SD * x).') / (obj.lamda + (SD * x).' * x ) ) / obj.lamda;
    w = w+SD*en(i)*x;
end
en = en(:);
yn = yn(:);
% 输出
ffe = yn;
end