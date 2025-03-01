function result = Comb(n, k)
% n 个不同元素中选择 k 个元素的组合数
%采用递归的思想
% Note:n>K
    if n == k
        result = 1;
    elseif k == 0
        result = 1;
    else
        result = Comb(n-1, k-1) + Comb(n-1, k);
    end
end