function result = ccdm_make_cfg_probability(d_sz, c_sz, p0)
    if nchoosek(c_sz, floor(c_sz/2)) <= (2^d_sz) 
           warning('c_sz=%d is too low to be used', c_sz);
    end
    one_bits_min = 1;
    while nchoosek(c_sz, one_bits_min) < 2^d_sz
        one_bits_min = one_bits_min + 1;
    end
    %选择 c_sz 个数中至少有一个一比特值的概率。
    p1 = 1-p0;
    %期望
    one_bits_prob = round(p1 * c_sz);
    %判断警告
    if one_bits_min > one_bits_prob
        warning('one_bits_min = %d > one_bits_prob = %d, falling to one_bits_min, p0=%.2f (desired p0=%.2f)', one_bits_min, one_bits_prob, 1-one_bits_min/c_sz, p0);
    end

    result = struct('one_bits', max(one_bits_min, one_bits_prob), 'd_sz', d_sz, 'c_sz', c_sz);
end