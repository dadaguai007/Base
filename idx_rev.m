function x = idx_rev(y, w, x_sz)
    x = zeros(1, x_sz);
    %将y转换为行向量
    if size(y,1)>size(y,2)
        y=y.';
    end
    %转换为十进制
    y_val = arr2val(y);
    x_sum = 0;
    for i = 1:x_sz
        n = x_sz - i;
        k = w - x_sum;
        
        if n < k
            C=0;
        else
%             C = Comb(n, k);
            C=nchoosek(n, k);
        end

        if y_val >= C
            x_sum = x_sum + 1;
            x(i) = 1;
            y_val = y_val - C;
        end
    end
end


function result = arr2val(arr)
%二进制转换为十进制
    result = sum(arr .* 2.^(0:(length(arr)-1)));
end