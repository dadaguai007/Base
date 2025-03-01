%% 创建对角矩阵
elements = [1, 2, 3, 4];
diagonal_matrix = diag(elements);
% 创建对角下一列矩阵
elements = [1, 2, 3];
lower_diag_matrix = diag(elements, -1);
sum_diag_martix=diagonal_matrix+lower_diag_matrix;
% 显示生成的对角矩阵
disp(diagonal_matrix);
% 打印生成的对角下一列矩阵
disp(lower_diag_matrix);

%打印生成求和矩阵
disp(sum_diag_martix);