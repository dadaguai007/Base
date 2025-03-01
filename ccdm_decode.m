function data_bits = ccdm_decode(code_bits, one_bits, d_sz, c_sz)
%     assert(mod(length(code_bits), c_sz) == 0, 'code_bits.size % c_sz != 0: %d %% %d = %d', length(code_bits), c_sz, mod(length(code_bits), c_sz));

    data_bits = zeros(1, floor(length(code_bits)/c_sz)*d_sz);
    
    for i = 1:floor(length(code_bits)/c_sz)
        data_bits((i-1)*d_sz+1:i*d_sz) = idx(code_bits((i-1)*c_sz+1:i*c_sz), one_bits, d_sz);
    end
end