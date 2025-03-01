function vnle2 = VNLE2_LMS(obj, xn, dn)
dn = dn(:);
L1 = length(xn);
L2 = length(dn);

% 输出长度
n = round(L1/obj.sps);
debugInfo = [];
xn = xn(:);

tap_ff = obj.k1;
tap_fb = obj.k2;

h1 = zeros(tap_ff,1);
x1 = zeros(tap_ff,1);
x2 = zeros(tap_fb,1);
% 参考抽头
x3 = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));
% x3 = xn(ref+1:end);

m = 1;
for k = 1:tap_fb
    for j = k:tap_fb
        index_mat(m,:) = [k,j];
        m = m+1;
    end
end
h2 = zeros(size(index_mat,1),1);

for i = 1:n - 2
    % 线性项
    x1 = cat(1,x1(obj.sps+1:end),x3(obj.sps*i-obj.sps+1:1:obj.sps*i));
    % 交叉项乘积
    x2 = x1(index_mat(:,1)).*x1(index_mat(:,2));
    yn(i,:) = x1.'*h1 + x2.'*h2;
    if i + tap_ff - 1 > L2
        %判决
        [~,y_d(i)] = quantiz(yn(i),[-2, 0, 2], [-3, -1, 1, 3]);
        en(i,:) = y_d(i)-yn(i);
    else
        en(i,:) = dn(i)-yn(i);
    end
    h1 = h1+obj.u*en(i)*x1;
    h2 = h2+obj.u*en(i)*x2;
    w = cat(1, h1, h2);
end
en = en(:);
yn = yn(:);
debugInfo.yn = yn;
vnle2 = yn;
end