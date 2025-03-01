function [dfe,en,w] = DFE_LMS(obj, xn, dn)
dn = dn(:);
L1 = length(xn);
L2 = length(dn);

% 抽头长度
k_fe = obj.k1;
k_fb = obj.k2;

n = round(L1/obj.sps);
yn = zeros(n, 1);
w1 = zeros(k_fe, 1);
w2 = zeros(k_fb, 1);
x1 = zeros(k_fe, 1);
x2 = zeros(k_fb, 1);
x3 = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));

for i = 1:n-1
    x1 = cat(1,x1(obj.sps+1:end),x3(obj.sps*i-obj.sps+1:1:obj.sps*i));
    yn(i) = x1.'*w1+ x2.'*w2;
    % 进入判决程序
    if i+ k_fe- 1 > L2
        if yn(i) > 2
            y_d(i) = 3;
        elseif yn(i) > 0
            y_d(i) = 1;
        elseif yn(i) > -2
            y_d(i) = -1;
        else
            y_d(i) = -3;
        end
        x2_new = y_d(i);
        en(i) = y_d(i)- yn(i);
    else
        x2_new = dn(i);
        en(i) = dn(i)- yn(i);
    end
    % 更新系数
    w1 = w1+ obj.u*en(i)*x1;
    w2 = w2+ obj.u*en(i)*x2;
    x2 = cat(1,x2(2:end),x2_new);
end
w = cat(1,w1,w2);
en = en(:);
yn = yn(:);

dfe = yn;
end