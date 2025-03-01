function vnle3full = VNLE3_full(obj, xn, dn)
L1 = length(xn);
L2 = length(dn);
% 输出长度
L = 2*floor(L1/obj.sps/2);

tap_nonlinear = floor(obj.tap/2);

x_linear = zeros(obj.tap, 1);
c1 = zeros(obj.tap, 1);
c2 = zeros(obj.tap*(obj.tap+1)/2, 1);
yn = zeros(L, 1);

% 参考抽头 % 按列进行串联
xn = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));

m = 1;
for i = 1:obj.tap
    for j = i:obj.tap
        index_mat_nonlinear1(m, :) = [i, j];
        m = m + 1;
    end
end


n = 1;
for ii = 1:tap_nonlinear
    for jj = ii:tap_nonlinear
        for kk = jj:tap_nonlinear
            index_mat_nonlinear2(n, :) = [ii, jj, kk];
            n = n+1;
        end
    end
end

c3 = zeros(size(index_mat_nonlinear2, 1), 1);

for idx = 1:L
    x_linear = cat(1,x_linear(obj.sps+1:end),xn(obj.sps*idx-obj.sps+1:1:obj.sps*idx));
    % 二阶交叉项
    x_nonlinear1 = x_linear(index_mat_nonlinear1(:, 1)).*x_linear(index_mat_nonlinear1(:, 2));
    % 三阶交叉项
    x_nonlinear2 = x_linear(index_mat_nonlinear2(:, 1)).*x_linear(index_mat_nonlinear2(:, 2)).*x_linear(index_mat_nonlinear2(:,3));
    % 均衡
    yn(idx) = x_linear.'*c1 + x_nonlinear1.'*c2 + x_nonlinear2.'*c3;
    if idx+ obj.tap- 1 > L2
        %判决函数
        [~,y_d(idx)] = quantiz(yn(idx),[-2, 0, 2], [-3, -1, 1, 3]);
        en(idx,:) = y_d(idx) - yn(idx);
    else
        en(idx,:) = dn(idx) - yn(idx);
    end
    c1 = c1 + obj.u1*en(idx)*x_linear;
    c2 = c2 + obj.u2*en(idx)*x_nonlinear1;
    c3 = c3 + obj.u3*en(idx)*x_nonlinear2;
end
yn = yn(:);
en = en(:);

vnle3full = yn;
end