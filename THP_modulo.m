function result = THP_modulo(x, A)
% 在DPC基础上进行取模的计算
    realPart = 2 * A * floor((real(x) + A) / (2 * A));
    imagPart = 2 * A * floor((imag(x) + A) / (2 * A));
    result = x - (realPart + 1i * imagPart);
end