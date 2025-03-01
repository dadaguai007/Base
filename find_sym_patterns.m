function [pattern_idx,sym_ptrns] = find_sym_patterns(sig, ref_sym, N, ret_ptrns)
% Find and index patterns elements of length N.
% LUT

if nargin < 4
    ret_ptrns = 'false'; % Default value for ret_ptrns if not specified
end
% 总共存在多少种符号
M = length(ref_sym);
%可能得模式数量
L = M^N;
%模式索引
idx=1:L;
% 将索引生成矩阵形式
idx = reshape(idx, (repmat(M, 1, N)));
% 将模式数组进行转置
% 循环数量为N-2
if N==3
    for i=1:(N+1)
        idx(:,:,i)=idx(:,:,i).';
    end
elseif N==4
    for j = 1:M
        for i=1:M
            idx(:,:,i,j)=idx(:,:,i,j).';
        end

    end
elseif N==5
    for jj= 1:M
        for j = 1:M
            for i=1:M
                idx(:,:,i,j,jj)=idx(:,:,i,j,jj).';
            end
        end
    end

end

% 广播减法生成矩阵
k=abs(bsxfun(@minus, sig(:).', ref_sym(:)));
%整个列表中在列中找最小值
[~, sig_idx]=min(k);

% 以N为一组，确定为一个模式 ,以发射信号为中心
sig_rwin = rolling_window_central(sig_idx, N, 'true');
sig_rwin=sig_rwin.';
% 第一行为哪一个矩阵，第二行为第几行，第三行为第几列(只适合N=3)
pattern_idx=zeros(length(sig_rwin),1);

% 选取相应的模式 一一对应
if N==3
    for i=1:length(sig_rwin)
        pattern_idx(i) = idx(sig_rwin(2,i),sig_rwin(3,i),sig_rwin(1,i));
    end

elseif N==4
    for i=1:length(sig_rwin)
        pattern_idx(i) = idx(sig_rwin(3,i),sig_rwin(4,i),sig_rwin(2,i),sig_rwin(1,i));
    end
elseif N==5
    for i=1:length(sig_rwin)
        pattern_idx(i) = idx(sig_rwin(4,i),sig_rwin(5,i),sig_rwin(3,i),sig_rwin(2,i),sig_rwin(1,i));
    end
end

if strcmp(ret_ptrns,'true')
    % maybe have the problem
    pidx = reshape(ind2sub(size(idx), 1:L), [], N);
    sym_ptrns = ref_sym(pidx);
else
    % 输出
    sym_ptrns = [];

end

end