function [rho, sigmasq] = yule_walker(X, order)
    % 计算 Yule-Walker 自回归系数。


    % 检查 X 是否为 2D 矩阵
%     assert(ndims(X) == 2, '输入矩阵必须为2D。');

    % 计算分母
    denom = size(X, 2) - (0:order);

    % 初始化自相关系数数组 r
    r = zeros(1, order + 1);

    % 遍历输入矩阵的每一行
    for di = 1:size(X, 1)
        % 从每一行中去除均值
        d = X(di, :) - mean(X(di, :));

        % 计算 lag 为 0 的自相关系数
        r(1) = r(1) + sum(d.^2);

        % 计算 lag 为 1 到 order 的自相关系数
        for k = 1:order
            r(k + 1) = r(k + 1) + sum(d(1:end - k) .* d(k + 1:end));
        end
    end

    % 将每个自相关系数标准化，即除以分母和样本长度
    r = r ./ (denom * size(X, 1));

    % 解 Yule-Walker 方程以得到自相关系数 rho
    rho = linsolve(toeplitz(r(1:end-1)), r(2:end)');

    % 计算噪声方差
    sigmasq = r(1) - sum(r(2:end) .* rho);

    % 取噪声方差的平方根
    sigmasq = sqrt(sigmasq);
end

% 调用该AR滤波器
% do the fitting
% coeffs = yule_walker(data, order);
% b = 1.0; % Numerator filter coefficients ， 分子系数
% a = [1.0, -coeffs]; % Denominator filter coefficients
%指递归滤波器（IIR filter）的分母系数。
% 应用滤波器
% y = filter(b, a, data);
% 创新序列  创新通常指的是模型预测值与实际观测值之间的差异。
% 在这里，创新是通过将信号 d 与递归滤波器的分母系数 a 进行卷积得到的。
% innovation = conv(d, a, "valid");
% 重新生成信号
% d_ = filter(b, a, innovation);
% 添加虚拟样本
% d_ = [ones(1, order)*d_(1), d_];
