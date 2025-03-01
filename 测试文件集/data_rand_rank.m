
%随机打乱函数的应用
% 一个已知的数据集 data，这里只是一个示例数据集
data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

% 使用 randperm 对数据集进行随机化
random_order = randperm(length(data));

% 使用随机顺序重新排列数据集
shuffled_data = data(random_order);

% 打印随机化后的数据集
disp(shuffled_data);