
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