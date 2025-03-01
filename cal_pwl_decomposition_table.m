function [M,debugInfo] = cal_pwl_decomposition_table(sig,lambda)
% strictly follow the table in the paper
% Y. Fu, JLT 2018, doi:10.1109/JLT.2019.2948096
debugInfo.IOI = []; 
sig = sig(:);
lambda = lambda(:);
% 阈值排序
lambda = sort(lambda);
N = length(lambda);
M = nan(length(sig),N+1);
%创建一个包含负无穷大、排序后的阈值向量和正无穷大的数组，用于构建一个查找表
tmp = [-inf;lambda;inf];
% 创建出N+1个区间，而且要包括一个正无穷和一个负无穷的区间阈值
% 所以，上面先加入两个无穷值，再使用buffer函数创建N+1个区间中每个区间的两边的阈值
% 输出数列应该为2*N+1 的样式
I = buffer(tmp,2,1,'nodelay'); % I has N+1 columns
for idx = 1:length(sig)
    gamma = sig(idx);
    for ii = 1:N+1
        lambda_i = I(2,ii);
        lambda_i_1 = I(1,ii);
        if gamma >= 0
            if gamma >= lambda_i_1 && gamma <= lambda_i
                if lambda_i_1 > 0
                    y = gamma - lambda_i_1;
                elseif lambda_i_1 <= 0 && lambda_i >= 0
                    y = gamma;
                end
            elseif gamma > lambda_i
                if lambda_i_1 > 0
                    y = lambda_i - lambda_i_1;
                    debugInfo.IOI(end+1) = idx;
                elseif lambda_i_1 <= 0 && lambda_i >= 0
                    y = lambda_i;
                elseif lambda_i < 0
                    y = 0;
                end
            elseif gamma <= lambda_i_1
                if lambda_i_1 > 0
                    y = 0;
                end
            end
        else % gamma < 0
            if gamma >= lambda_i_1 && gamma <= lambda_i
                if lambda_i_1 <= 0 && lambda_i >= 0
                    y = gamma;
                else
                    y = gamma - lambda_i;
                end
            elseif gamma > lambda_i
                if lambda_i < 0
                    y = 0;
                end
            elseif gamma <= lambda_i_1
                if lambda_i_1 > 0
                    y = 0;
                elseif lambda_i_1 <=0 && lambda_i >= 0
                    y = lambda_i_1;
                elseif lambda_i < 0
                    y = lambda_i_1 - lambda_i;
                end
            end
        end
       M(idx,ii) = y; 
    end
end



%
% 这段代码使用MATLAB的`buffer`函数来创建一个矩阵 `I`，其中存储了输入向量 `tmp` 中相邻两个元素的组合。让我们详细解释这行代码：
% I = buffer(tmp, 2, 1, 'nodelay'); 
% - `tmp`: 输入的向量，这里是`[-inf; lambda; inf]`。
% - `2`: 表示每一列包含的元素数目。在这里，每一列包含两个相邻的元素。
% - `1`: 表示相邻列之间的偏移量。在这里，相邻列之间没有重叠。
% - `'nodelay'`: 表示不延迟数据。这意味着第一列从向量的第一个元素开始，而不是从第二个元素开始。
% 
% 因此，这行代码的目的是将输入向量 `tmp` 中的相邻两个元素组合成一个矩阵 `I`，其中每一列包含两个相邻的元素。在这个特定的上下文中，`I` 用于存储阈值的相邻对，方便后续的计算。