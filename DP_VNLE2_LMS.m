% DP-VNLE  这些VNLE都是FFE形状的，并没有进行反馈，也就不是DFE结构
function pnle2 = DP_VNLE2_LMS(obj, xn, dn)

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

tap_fb = obj.k2;  


% 生成上三角矩阵的系数矩阵
upper_triangular_matrix = zeros(tap_fb);
% 选取对角线条数
W=obj.W;
% 根据对角线的数量，填充对角线元素
for i = 1:n
    for j=i:i+W-1
        if j>n
            break
        end
    upper_triangular_matrix(i, j) = j;
    end
end


for i = 1:L
    x_linear = cat(1,x_linear(obj.sps+1:end),xn1(obj.sps*i-obj.sps+1:1:obj.sps*i));
    % 选取W条对角线交叉项 注意sum求和，两个矩阵相乘
    x_nonlinear = sum(x_linear(1:obj.k2).*upper_triangular_matrix,2);
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