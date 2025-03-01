function vnle3 = VNLE3_LMS(obj, xn, dn)
dn = dn(:);
L1 = length(xn);
L2 = length(dn);
debugInfo = [];

% 信号输出长度
n = round(L1/obj.sps);
xn = xn(:);

yn = zeros(n, 1);
w1 = zeros(obj.k1, 1);
x1 = zeros(obj.k1, 1);
x2 = zeros(obj.k2, 1);
x3 = zeros(obj.k3, 1);
% 参考抽头
x4 = cat(1, xn(obj.ref+1:end), dn(end-obj.ref+1:end));

% 生成volterra矩阵
m = 1;
for k = 1:obj.k2
    for j = k:obj.k2
        for q = j:obj.k3
            index_mat(m,:) = [k,j,q];
            m = m+1;
        end
    end
end
w2 = zeros(size(index_mat,1),1);
w3 = zeros(size(index_mat,1),1);

for i = 1:n - 1
    x1 = cat(1,x1(obj.sps+1:end),x4(obj.sps*i-obj.sps+1:1:obj.sps*i));
    %二阶交叉项
    x2 = x1(index_mat(:,1)).*x1(index_mat(:,2));
    % 三阶交叉项
    x3 = x1(index_mat(:,1)).*x1(index_mat(:,2)).*x1(index_mat(:,3));

    yn(i,:) = x1.'*w1 + x2.'*w2 + x3.'*w3;
    if i > L2
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
    w1 = w1+obj.u1*en(i)*x1;
    w2 = w2+obj.u2*en(i)*x2;
    w3 = w3+obj.u3*en(i)*x3;
    w = cat(1, w1, w2, w3);
end
en = en(:);
yn = yn(:);
vnle3 = yn;
end