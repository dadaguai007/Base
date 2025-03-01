function iq = qam_mapping(b, M)
% 自己写的qam调制代码
b=b.';
bits = log2(M);
    if mod(length(b), bits) ~= 0
        warning('b.size not divisible by %d, zeros appended', bits);
        b = [b, zeros(1, bits - mod(length(b), bits))];
    end

    % get chosen constellation map
    sym2iq_map = qammod((0:M-1)',M,'gray');

    % generate iq points based on chosen constellation map
    idx = zeros(1, length(b)/bits);
    for i = 1:bits
        %二进制左移
        idx = idx + bitshift(b(i:bits:end), i-1);
    end

    iq = sym2iq_map(idx + 1);
end