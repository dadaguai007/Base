function wout = my_nl_mlse(win, M, const, ML, lcoeff, nlcoeff)
% 示例：% mlse_out = my_nl_mlse(ffe_out, M, constellation, taps1, lcoeff, []).';
% alpha_pnc = 0.4;
%pf_out = mlseeq(filter([1, alpha_pnc], 1, ffe_out), [1, alpha_pnc], constellation, 1000, 'rst');
if nargin < 6
    nlcoeff = [];
elseif nargin ~= 6
    error('wout = my_nl_mlse(win, M, const, lcoeff, nlcoeff, ML)');
end
[lut, ref] = genlut(M, const, ML, lcoeff, nlcoeff);
win = win(:).';
wout = zeros(1, length(win));
pm = zeros(1, M^(ML(1)));% 状态量
for idx = 1:length(win)
    pmVec = pm - (win(idx)-lut).*conj(win(idx)-lut);
    [~, sind] = max(pmVec);  
    pmVec = reshape(pmVec, M, []);
    wout(idx) = ref(1, sind); % 找到接近原始向量的值，sind代表置信度最大的路径

    [pm,~] = max(pmVec);
    pm = repmat(pm,1,M); % refer to Sheet2@Excel file
end
% 除去开头处的溯源态
wout(1:ML(1)-1) = [];
end

% 实现了一种基于查找表（LUT）的非线性最大似然序列估计（MLSE）算法。
% ```
% - **输入参数**：
%     - `win`：输入信号向量。
%     - `M`：调制阶数，例如对于QPSK，`M = 4`。
%     - `const`：星座图，用于定义调制符号的位置。
%     - `ML`：记忆长度等。
%     - `lcoeff`：线性滤波器系数。
%     - `nlcoeff`：非线性滤波器系数（如果未提供，将初始化为空数组）。
% - **输出参数**：
%     - `wout`：经过非线性最大似然序列估计处理后的输出信号向量。
% 
% ### 参数检查与初始化
% - 如果输入参数个数小于6，将`nlcoeff`初始化为空数组。
% - 如果输入参数个数不等于6，抛出错误并提示正确的函数调用格式。
% 
% ### 生成查找表和参考向量
% ```matlab
% [lut, ref] = genlut(M, const, ML, lcoeff, nlcoeff);
% ```
% - 调用`genlut`函数（该函数需自定义实现），
% 根据输入的参数`M`、`const`、`ML`、`lcoeff`和`nlcoeff`生成查找表`lut`和参考向量`ref`。
% 这个查找表和参考向量是后续MLSE算法的关键数据。

% ### 信号预处理和初始化
% ```matlab
% win = win(:).';
% wout = zeros(1, length(win));
% pm = zeros(1, M^(ML(1)));
% ```
% - `win = win(:).';`：将输入信号向量`win`转换为行向量。
% - `wout = zeros(1, length(win));`：初始化输出信号向量`wout`为全零向量，长度与输入信号`win`相同。
% - `pm = zeros(1, M^(ML(1)));`：初始化一个概率度量向量`pm`，长度为`M`的`ML(1)`次幂，用于存储每个可能状态的概率度量。
% 
% ### 非线性最大似然序列估计循环

% - 遍历输入信号向量`win`的每个元素：
%     - `pmVec = pm - (win(idx)-lut).*conj(win(idx)-lut);`：计算当前输入信号与查找表中每个元素的误差平方和，
% 并从当前的概率度量向量`pm`中减去，得到新的概率度量向量`pmVec`。
%     - `[~, sind] = max(pmVec);`：找到`pmVec`中最大值的索引`sind`。
%     - `pmVec = reshape(pmVec, M, []);`：将`pmVec`重塑为一个`M`行的矩阵，以便后续处理。
%     - `wout(idx) = ref(1, sind);`：根据索引`sind`从参考向量`ref`中获取对应的估计值，并赋值给输出信号向量`wout`的当前位置。
%     - `[pm,~] = max(pmVec);`：找到重塑后的`pmVec`中每一列的最大值，更新概率度量向量`pm`。
%     - `pm = repmat(pm,1,M);`：将更新后的`pm`向量重复`M`次，这一步可能是为了与后续的计算或算法逻辑相匹配，

% 
% ### 输出信号处理
% - 由于MLSE算法可能存在初始的暂态过程，这里将输出信号向量`wout`的前`ML(1) - 1`个元素删除，得到最终的输出信号。
% 
% 
% 总体而言，这段代码实现了一个基于查找表的非线性最大似然序列估计算法，通过对输入信号的逐元素处理，利用查找表和概率度量更新来估计输出信号。