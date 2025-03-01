%% 使用概率分布来实现星座整形
% 理想情况下的实现


function [p0, p1] = qam_optimal_p0(C)
%     alpha = fminbnd(@(x) qam_C_optim(x, C), 0, 1);
    % 定义待优化的函数
    fun = @(x) qam_C_optim(x, C);

    % 设置初始猜测值的区间
    interval = [0, 1];

    % 使用 fzero 寻找根
    alpha = fzero(fun, interval);
    % 生成比特0的概率
    p0 = exp(-alpha);
    %生成比特1的概率
    p1 = exp(-9*alpha);
    norm = p0 + p1;
    p0 = p0 / norm;
    p1 = p1 / norm;
end

function res = qam_C_optim(alpha, C)
    pts_sqr = [2, 10, 10, 18];
    p_pts = exp(-alpha * pts_sqr);
    p_pts = p_pts ./ sum(p_pts);

    res = -C - sum(p_pts .* log2(p_pts)) + 2;
end