function [Error, idx_I, idx_Q] = cal_lut(tx_sig, rx_sig, ref_sym, mem_len, idx_data, real_ptrns)
% Calculate a lookup table for precompensation of pattern-based errors.
%基于模式的误差预补偿的查找表。可以视为简化的 Volterra 过滤器。仅适用于 1D 信号，数据有更多的维度，需要循环使用。
% mem_len 补偿长度
% real_ptrns 是否使用复数数组变换
%如果为实质处理，输出的两项相等
    if nargin < 5
        idx_data = ones(size(tx_sig));
    end
    if nargin < 6
        real_ptrns = 'true';
    end
    err = (tx_sig - rx_sig).';

    idx = find(idx_data>0); 


    if strcmp(real_ptrns,'false')
        %去掉重复值
        ref_sym_I = unique(real(ref_sym));
        ref_sym_Q = unique(imag(ref_sym));
        M = length(ref_sym_I);
        N = M.^mem_len;

        idx_I = find_sym_patterns(real(tx_sig), ref_sym_I, mem_len);
        idx_Q = find_sym_patterns(imag(tx_sig), ref_sym_Q, mem_len);
        idx_I=idx_I(idx);
        idx_Q=idx_Q(idx);
        % 建立LUT
        ea = cal_lut_avg(err, idx_I, idx_Q, N);
        real_part=real(ea);
        imag_part=imag(ea);
        Error=real_part(idx_I).'+1j*imag_part(idx_I).';
        idx_I = idx_I.';
        idx_Q = idx_Q.';
    elseif strcmp(real_ptrns,'true')
        %去掉重复值
        ref_sym_c = unique(ref_sym);
        M = length(ref_sym_c);
        N = M.^mem_len;
        % 找到模式
        idx_c = find_sym_patterns(tx_sig, ref_sym_c, mem_len);
        idx_c=idx_c(idx);
        % 建立LUT
        ea = cal_lut_avg(err, idx_c, idx_c, N);
        Error=ea(idx_c).';
        idx_I = idx_c.';
        idx_Q = idx_c.';
    end


end