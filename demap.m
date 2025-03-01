
function decBits = demap(indSymb, bitMap)
% Contellation symbol index to bit sequence demapping.
%
% Parameters
% ----------
% indSymb : np.array of ints
%     Indexes of received symbol sequence.
% bitMap : (M, log2(M)) np.array
%     bit-to-symbol mapping.
%
% Returns
% -------
% decBits : np.array
%     Sequence of demapped bits.

M = size(bitMap, 1);
b = log2(M);

decBits = zeros(length(indSymb) * b, 1);

for i = 1:length(indSymb)
    decBits(i * b - b + 1:i * b) = bitMap(indSymb(i) + 1, :); % Matlab index starts from 1
end
end

%这段代码是一个名为`demap`的函数，用于将星座符号的索引序列映射回对应的比特序列。这是无线通信中调制解调过程的一部分，特别是在数字通信系统中。函数使用两个参数：`indSymb`和`bitMap`。
% 参数解释如下：
% - `indSymb`：这是一个整数类型的NumPy数组，表示接收到的符号序列的索引。每个索引对应于星座图中的一个点，每个点代表一个或多个比特。
% - `bitMap`：这是一个二维NumPy数组，其大小为`M x log2(M)`，其中`M`是星座点数。这个数组定义了从比特到符号的映射。每一行代表一个符号，对应的比特序列按照列的顺序排列。
% 函数的返回值是`decBits`，这是一个NumPy数组，包含了解映射后的比特序列。
% 函数的工作流程如下：
% 1. `M = size(bitMap, 1);`获取`bitMap`数组的行数，即星座点数`M`。
% 2. `b = log2(M);`计算每个符号对应的比特数，这是通过对星座点数取以2为底的对数得到的。
% 3. `decBits = zeros(length(indSymb) * b, 1);`初始化解映射后的比特序列数组，其长度为`indSymb`数组长度乘以每个符号的比特数。
% 4. 循环遍历`indSymb`数组中的每个索引：
%    - `i = 1:length(indSymb)`是循环的迭代变量。
%    
%    - `decBits(i * b - b + 1:i * b) = bitMap(indSymb(i) + 1, :);`这行代码将`indSymb`中的索引映射到`bitMap`中对应的比特序列，并将其存储在`decBits`数组的正确位置。注意，由于MATLAB的索引是从1开始的，所以这里`indSymb(i) + 1`是为了匹配MATLAB的索引习惯。
% 最终，`decBits`数组包含了原始比特序列，可以用于后续的解码和处理。
