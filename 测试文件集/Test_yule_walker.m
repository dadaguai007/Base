% 假设我们有一个时间序列数据集X，这里我们使用随机数生成一个模拟的AR(2)过程作为示例
order = 2; % 我们模拟一个二阶自回归过程
n = 100; % 时间序列的长度
phi1 = -0.3; % 第一个自回归系数
phi2 = 0.7; % 第二个自回归系数
sigma = 0.5; % 噪声项的标准差

% 生成模拟数据
e = sigma * randn(n, 1); % 噪声项
X = zeros(n, 1);
X(1) = e(1); % 初始化第一个值
X(2) = phi1 * X(1) + phi2 * X(1) + e(2); % 初始化第二个值

for i = 3:n
    X(i) = phi1 * X(i - 1) + phi2 * X(i - 2) + e(i);
end

% 测试yule_walker函数
[rho, sigmasq] = yule_walker(X, order);

% 显示结果
disp('估计的自回归系数 rho:');
disp(rho);
disp('估计的噪声方差 sigmasq:');
disp(sigmasq);

% 可选：可视化结果
figure;
subplot(2,1,1);
plot(X);
title('Simulated AR(2) Process');
subplot(2,1,2);
plot(rho);
title('Estimated AR Coefficients');