function result = idx(x, w, y_sz, warnings_on)
    if nargin < 4
        warnings_on = false;
    end

    if w ~= sum(x)
        if warnings_on
            warning("Bad Sequence w != sum(x): %d != %d, ignoring %d MSB ones", w, sum(x), sum(x)-w);
        end

        to_ignore = sum(x) - w;
        for i = 1:length(x)
            if x(end-i+1) == 1
                x(end-i+1) = 0;
                to_ignore = to_ignore - 1;
                if to_ignore == 0
                    break;
                end
            end
        end
    end

    res = 0;
    x_sum = 0;
    for i = 1:length(x)
        if x(i) ~= 0
            n = length(x) - i;
            k = w - x_sum;
            x_sum = x_sum + x(i);
            if k <= n
                res = res + nchoosek(n, k);
            end
        end
    end
    result = val2arr(res, y_sz);
end


function [res,Res] = val2arr(val, bits)
% 将整数转换为二进制，当然转换后的二进制结果需要进行一个前后颠倒
% 第一位位于索引1处
    res = zeros(1, bits);
    for i = 1:bits
        res(i) = bitand(bitshift(val, -(i-1)), 1);
    end
    Res=flip(res);
end