function res = val2arr(val, bits)
    res = zeros(1, bits);
    for i = 1:bits
        res(i) = bitand(bitshift(val, -(i-1)), 1);
    end
end

function result = arr2val(arr)
    result = sum(arr .* 2.^(0:(length(arr)-1)));
end





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
                res = res + Comb(n, k);
            end
        end
    end
    result = val2arr(res, y_sz);
end

function x = idx_rev(y, w, x_sz)
    x = zeros(1, x_sz);

    y_val = arr2val(y);
    x_sum = 0;
    for i = 1:x_sz
        n = x_sz - i;
        k = w - x_sum;
        C = 0;
        if n >= k
            C = Comb(n, k);
        end

        if y_val >= C
            x_sum = x_sum + 1;
            x(i) = 1;
            y_val = y_val - C;
        end
    end
end


function code_bits = ccdm_encode(data_bits, one_bits, d_sz, c_sz)
    
    if mod(length(data_bits), d_sz) ~= 0
        warning('data.size not divisible by d_sz: %d%%%d=%d, zeros appended', length(data_bits), d_sz, mod(length(data_bits), d_sz));
        data_bits = [data_bits, zeros(1, mod(length(data_bits), d_sz))];
    end

    code_bits = zeros(1, floor(length(data_bits)/d_sz)*c_sz);

    for i = 1:floor(length(data_bits)/d_sz)

        code_bits((i-1)*c_sz+1:i*c_sz) = idx_rev(data_bits((i-1)*d_sz+1:i*d_sz), one_bits, c_sz);
    
    end
end

function data_bits = ccdm_decode(code_bits, one_bits, d_sz, c_sz)
    assert(mod(length(code_bits), c_sz) == 0, 'code_bits.size % c_sz != 0: %d %% %d = %d', length(code_bits), c_sz, mod(length(code_bits), c_sz));

    data_bits = zeros(1, length(code_bits)/c_sz*d_sz);
    
    for i = 1:length(code_bits)/c_sz
        data_bits((i-1)*d_sz+1:i*d_sz) = idx(code_bits((i-1)*c_sz+1:i*c_sz), one_bits, d_sz);
    end
end








