function [lut, ref] = genlut(M, const, ML, lcoeff, nlcoeff)
% ML为回溯的记忆长度
if nargin < 5
    nlcoeff = [];
elseif nargin ~= 5
    error('genlut(M, const, lcoeff, nlcoeff, ML)');
end
if max(ML) ~= ML(1)
    error('The memory length of the linear section should be maximum.');
end
% 生成一个`L`行、`M^L`列的矩阵，
% 每一列代表一个长度为"回溯长度"的符号序列,生成状态矩阵的序号
idx = getSymVecFromML(M, ML(1));

% 根据序号生成状态矩阵
ref = const(idx);

w = [lcoeff, nlcoeff];
if isempty(nlcoeff)
    lut = w * ref;% 生成线性LUT，（线性向量与状态向量的乘积）
else
    block = ref;
    for j = 2:length(ML)
        tempIdx = ML(1)-ML(j)+1:ML(1);
        block = [block; gen_vol_block(block(tempIdx, :), j)];
    end
    lut = w * block;
end
end

% 
% 主要功能是生成查找表（LUT, Look-Up Table）和参考信号
% 
% - 该函数接受五个输入参数：
%     - `M`：可能是调制阶数或与信号相关的维度参数。
%     - `const`：星座点集合，代表了信号可能的取值。
%     - `ML`：记忆长度向量，用于定义线性和非线性部分的记忆长度。
%     - `lcoeff`：线性系数向量。
%     - `nlcoeff`：非线性系数向量（可选参数）。
% - 函数返回两个输出：
%     - `lut`：生成的查找表。
%     - `ref`：参考信号。
% 

% 
% ### 记忆长度检查
% ```matlab
% if max(ML) ~= ML(1)
%     error('The memory length of the linear section should be maximum.');
% end
% ```
% - 检查`ML`向量中的最大值是否等于第一个元素。如果不相等，抛出错误，提示线性部分的记忆长度应该是最大的。
% 
% ### 生成参考信号
% ```matlab
% idx = getSymVecFromML(M, ML(1));
% ref = const(idx);
% ```
% - 调用`getSymVecFromML`函数（该函数在当前代码中未定义），传入`M`和`ML(1)`，得到一个索引向量`idx`。
% - 使用索引向量`idx`从星座点集合`const`中提取相应的元素，生成参考信号`ref`。
% 
% ### 生成查找表
% ```
% - 将线性系数向量`lcoeff`和非线性系数向量`nlcoeff`合并成一个向量`w`。
% - 如果`nlcoeff`为空向量：
%     - 直接将`w`与`ref`相乘，得到查找表`lut`。
% - 如果`nlcoeff`不为空向量：
%     - 初始化`block`为`ref`。
%     - 通过循环遍历`ML`向量从第二个元素开始：
%         - 计算一个临时索引向量`tempIdx`。
%         - 调用`gen_vol_block`函数（该函数在当前代码中未定义），传入`block`的一部分和当前的索引`j`，生成一个新的块。
%         - 将新生成的块与`block`垂直拼接。
%     - 最后，将`w`与拼接后的`block`相乘，得到查找表`lut`。
% 