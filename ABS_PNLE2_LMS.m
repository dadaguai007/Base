% ABS-PNLE

function pnle2 = ABS_PNLE2_LMS(obj, xn, dn)

L1 = length(xn);
L2 = length(dn);
% 输出长度
L = 2*floor(L1/obj.sps/2);

x_linear = zeros(obj.k1, 1);
c1 = zeros(obj.k1, 1);
c2 = zeros(obj.k2, 1);
yn = zeros(L, 1);

% 参考抽头
xn1 = cat(1, xn(obj.ref+1:end), dn(end-obj.ref+1:end));

for i = 1:L
    x_linear = cat(1,x_linear(obj.sps+1:end),xn1(obj.sps*i-obj.sps+1:1:obj.sps*i));
    % 将非线性项改变为绝对值
    x_nonlinear = abs(x_linear(1:obj.k2));
    yn(i,:) = x_linear.'*c1 + x_nonlinear.'*c2;

    if i + obj.k1 -1 > L2
        % 判决
        if yn(i) > 2
            y_d(i) = 3;
        elseif yn(i) > 0
            y_d(i) = 1;
        elseif yn(i) > -2
            y_d(i) = -1;
        else
            y_d(i) = -3;
        end
        en(i,:) = y_d(i)-yn(i);
    else
        en(i,:) = dn(i)-yn(i);
    end
    % LMS更新抽头
    c1 = c1 + en(i)*obj.u1*x_linear;
    c2 = c2 + en(i)*obj.u2*x_nonlinear;
end
yn = yn(:);
en = en(:);
pnle2 = yn;
end