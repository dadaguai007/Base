function ref=get_state(M, const, ML)
% 生成序列的状态向量，用于构建状态和分支的符号序列

if max(ML) ~= ML(1)
    error('The memory length of the linear section should be maximum.');
end
% 生成一个`L`行、`M^L`列的矩阵，
% 每一列代表一个长度为"回溯长度"的符号序列,生成状态矩阵的序号
idx = getSymVecFromML(M, ML(1));

% 根据序号生成状态矩阵
ref = const(idx);

end