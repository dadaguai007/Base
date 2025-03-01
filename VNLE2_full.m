function vnle2_full = VNLE2_full(obj, xn, dn)

L1 = length(xn);
L2 = length(dn);
% 输出长度
L = 2* floor(L1/obj.sps/2);

x_linear = zeros(obj.tap, 1);
% 抽头
w1 = zeros(obj.tap, 1);
w2 = zeros(obj.tap*(obj.tap+1)/2, 1);
yn = zeros(L, 1);

% 参考抽头
xn = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));
% 设置好相应的阶乘项
m = 1;
for i = 1:obj.tap
    for j = i:obj.tap
        index_mat(m, :) = [i, j];
        m = m+1;
    end
end

for ii = 1:L
    x_linear = cat(1,x_linear(obj.sps+1:end),xn(obj.sps*ii-obj.sps+1:1:obj.sps*ii));
    x_nonlinear =  x_linear(index_mat(:, 1)).*x_linear(index_mat(:, 2));
    yn(ii) = x_linear.'*w1 + x_nonlinear.'*w2;
    if ii+ obj.tap- 1> L2
        % 判决
        [~,y_d(ii)] = quantiz(yn(ii),[-2, 0, 2], [-3, -1, 1, 3]);
        en(ii,:) = y_d(ii) - yn(ii);
    else
        en(ii,:) = dn(ii) - yn(ii);
    end
    w1 = w1+ 2*obj.u1* en(ii)* x_linear;
    w2 = w2+ 2*obj.u2* en(ii)* x_nonlinear;
end

yn = yn(:);
en = en(:);

vnle2_full = yn;
end