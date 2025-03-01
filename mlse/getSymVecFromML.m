function retIndMat = getSymVecFromML(M,L)
% this function generate L-by-M^L matrix that list all the possible
% combination of [I1,I2,...,IL] where I belongs to set [1:M]
% example:% created on 2023/07/12 by Tianwai@OCG,BIT
%   M = 2; L = 3
%   expected output:
%       111,112,121,122,211,212,221,222
% created on 2023/07/12 by Tianwai@OCG,BIT
mat = repmat(1:M,L,1);
c = mat2cell(mat,ones(1,L));
[x{1 : numel(c)}] = ndgrid(c{:});
retIndMat = reshape(cat(numel(c), x{:}), [], numel(c));
retIndMat = retIndMat.';
end


% 主要功能是生成一个包含所有可能符号组合的矩阵。
% 具体来说，生成一个`L`行、`M^L`列的矩阵，其中每一列代表一个长度为`L`的符号序列，每个符号都来自集合`[1:M]`。
% **输入参数**
%     - `M`：一个正整数，表示符号集合的大小。例如，当`M = 2`时，符号集合为`{1, 2}`。
%     - `L`：一个正整数，表示符号序列的长度。
% **函数主体详细分析**
%     - **创建基础矩阵**：
%       mat = repmat(1:M,L,1);

%       使用`repmat`函数创建一个`L`行、1列的矩阵`mat`，其中每一行都是从1到`M`的序列。`repmat`函数用于复制和排列矩阵，这里将向量`1:M`复制`L`次，形成一个`L`行的矩阵。

%     - **将矩阵转换为元胞数组**：
%       ```matlab
%       c = mat2cell(mat,ones(1,L));
%       ```
%       使用`mat2cell`函数将矩阵`mat`转换为元胞数组`c`。
% `mat2cell`函数根据指定的行和列划分方式将矩阵转换为元胞数组，
% 这里`ones(1,L)`表示将`mat`的每一行作为一个单独的元胞元素。

%     - **生成所有可能的组合**：
%       ```matlab
%       [x{1 : numel(c)}] = ndgrid(c{:});
%       ```
%       使用`ndgrid`函数生成所有可能的组合。`ndgrid`函数接受多个输入数组，
% 并生成这些数组的所有可能组合的网格。
% 输入的是元胞数组`c`中的各个元素，`numel(c)`表示元胞数组`c`的元素个数，即`L`。
% 生成的`x`是一个包含`L`个元素的元胞数组，每个元素都是一个包含所有可能组合的矩阵。

%     - **重塑和转置结果**：
%       ```matlab
%       retIndMat = reshape(cat(numel(c), x{:}), [], numel(c));
%       retIndMat = retIndMat.';
%       ```
%       首先，使用`cat`函数将`x`中的所有元胞数组按维度`numel(c)`（即`L`）连接起来，形成一个大的矩阵。
% 
% 然后，使用`reshape`函数将这个矩阵重塑为一个列数为`numel(c)`（即`L`），行数为`M^L`的矩阵。最后，对结果矩阵进行转置，得到最终的输出矩阵`retIndMat`，其中每一列代表一个长度为`L`的符号序列。

%     - 函数返回一个`L`行、`M^L`列的矩阵`retIndMat`，其中包含了所有可能的符号组合。
% 
% 例如，当`M = 2`，`L = 3`时：
% - `mat = repmat(1:2,3,1)`生成：

% - `c = mat2cell(mat,ones(1,3))`将其转换为元胞数组`c`，`c`包含三个元胞，每个元胞是一个`1x2`的向量`[1 2]`。
% - `[x{1 : 3}] = ndgrid(c{:})`生成所有可能的组合，`x`中的每个元素是一个矩阵：
% ```
% - `retIndMat = reshape(cat(3, x{:}), [], 3)`将这些矩阵连接并重塑为：
% ```
% 1 1 1 1 2 2 2 2
% 1 1 2 2 1 1 2 2
% 1 2 1 2 1 2 1 2
% ```

